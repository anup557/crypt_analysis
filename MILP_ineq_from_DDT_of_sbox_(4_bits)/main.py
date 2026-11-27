"""
MILP Linear Inequality Generation from Difference Distribution Table (DDT)
===========================================================================
Pipeline:
    1. Input S-box lookup table.
    2. Compute the Difference Distribution Table (DDT) of the S-box.
    3. Extract all valid differential transitions (DDT[in_diff][out_diff] > 0)
       as 2n-dimensional binary vectors (points in {0, 1}^2n).
    4. Compute H-Representation (Convex Hull facet hyperplanes) in exact
       coprime integer coefficients (A * x + b >= 0) using SciPy ConvexHull.
    5. Reduce inequalities using Algorithm 1 (Greedy Elimination) to obtain
       the minimal subset of linear inequalities modeling the valid DDT transitions.
    6. Verify that 100% of valid transitions satisfy all reduced inequalities,
       and all invalid transitions violate at least one inequality.
"""

import os
import sys
import math

from convex_hull_hrep import compute_h_representation, format_algebraic
from reduced import reduce_inequalities


def find_ddt(sbox):
    """
    Computes the Difference Distribution Table (DDT) of an n-bit S-box.
    Returns:
        ddt: 2^n x 2^n integer matrix of differential transition counts.
        n: S-box bit size.
    """
    size = len(sbox)
    n = int(math.log2(size))
    assert 2**n == size, f"S-box length ({size}) must be a power of 2!"

    ddt = [[0 for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for in_diff in range(size):
            out_diff = sbox[row] ^ sbox[row ^ in_diff]
            ddt[in_diff][out_diff] += 1

    return ddt, n


def print_ddt(ddt, cipher_name="S-box"):
    """
    Prints the DDT in formatted grid with '.' for zero entries.
    """
    size = len(ddt)
    print("=" * 80)
    print(f"DIFFERENCE DISTRIBUTION TABLE (DDT) FOR {cipher_name.upper()} ({size}x{size})")
    print("=" * 80)
    print("i\\o |  ", end="")
    idx_list = [f"{i:x}" for i in range(size)]
    for j in range(size):
        print(f"{idx_list[j]:>2}", end=" ")
    print("\n" + "-" * (7 + 3 * size))

    for i in range(size):
        print(f"{idx_list[i]:>2}  |  ", end="")
        for j in range(size):
            if ddt[i][j] == 0:
                print(" .", end=" ")
            else:
                print(f"{ddt[i][j]:>2}", end=" ")
        print()
    print()


def extract_valid_transitions(ddt, n):
    """
    Extracts all valid differential transitions (in_diff -> out_diff) where DDT > 0.
    Each transition is converted to a 2n-bit binary vector:
        [x_{n-1}, ..., x_0, y_{n-1}, ..., y_0]
    where x is in_diff (MSB to LSB) and y is out_diff (MSB to LSB).
    """
    size = 2**n
    valid_transitions = []
    for in_diff in range(size):
        for out_diff in range(size):
            if ddt[in_diff][out_diff] > 0:
                in_bits = [(in_diff >> i) & 1 for i in reversed(range(n))]
                out_bits = [(out_diff >> i) & 1 for i in reversed(range(n))]
                valid_transitions.append(in_bits + out_bits)
    return valid_transitions


def verify_inequalities(valid_points, reduced_ineqs, n):
    """
    Rigorously verifies:
      1. Every valid DDT point satisfies all reduced inequalities (>= 0).
      2. Every invalid DDT point violates at least one reduced inequality (< 0).
    """
    valid_set = set(tuple(p) for p in valid_points)
    total_binary_points = 2 ** (2 * n)
    all_points = [[(i >> k) & 1 for k in reversed(range(2 * n))] for i in range(total_binary_points)]

    # Check 1: Validity of valid points
    for p in valid_points:
        for ineq in reduced_ineqs:
            val = sum(x * y for x, y in zip(p, ineq[:-1])) + ineq[-1]
            if val < 0:
                return False, f"False Negative: Valid transition {p} violated inequality {ineq} (val={val})"

    # Check 2: Rejection of invalid points
    for p in all_points:
        if tuple(p) not in valid_set:
            violated = any((sum(x * y for x, y in zip(p, ineq[:-1])) + ineq[-1]) < 0 for ineq in reduced_ineqs)
            if not violated:
                return False, f"False Positive: Invalid transition {p} was NOT rejected by any inequality!"

    return True, f"All {len(valid_points)} valid transitions satisfied, and all {total_binary_points - len(valid_points)} invalid transitions eliminated!"


def main():
    # --------------------------------------------------------------------------
    # S-box Selection
    # --------------------------------------------------------------------------
    ciphers = {
        "DEFAULT": [0x0, 0x3, 0x7, 0xe, 0xd, 0x4, 0xa, 0x9, 0xc, 0xf, 0x1, 0x8, 0xb, 0x2, 0x6, 0x5],
        "PRESENT": [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2],
        "GIFT":    [0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe],
        "BAKSHEESH": [0x3, 0x0, 0x6, 0xD, 0xB, 0x5, 0x8, 0xE, 0xC, 0xF, 0x9, 0x2, 0x4, 0xA, 0x7, 0x1],
        "SLIM":    [0xc, 0x5, 0x6, 0xb, 0x9, 0x0, 0xa, 0xd, 0x3, 0xe, 0xf, 0x8, 0x4, 0x7, 0x1, 0x2],
        "FUTURE":  [0x1, 0x3, 0x0, 0x2, 0x7, 0xe, 0x4, 0xd, 0x9, 0xa, 0xc, 0x6, 0xf, 0x5, 0x8, 0xb],
        "EQ_292":  [0x8, 0x0, 0xA, 0x2, 0x3, 0x5, 0x1, 0x7, 0x4, 0x9, 0x6, 0xB, 0xC, 0xD, 0xE, 0xF]
    }

    # Choose cipher to model (can pass via command line argument, e.g. python3 ddt_main.py PRESENT)
    cipher_name = "DEFAULT"
    if len(sys.argv) > 1 and sys.argv[1].upper() in ciphers:
        cipher_name = sys.argv[1].upper()
    elif len(sys.argv) > 1:
        print(f"Unknown cipher '{sys.argv[1]}'. Available: {list(ciphers.keys())}")
        sys.exit(1)

    sbox = ciphers[cipher_name]

    print("=" * 80)
    print(f"DDT LINEAR INEQUALITY MODELING PIPELINE FOR {cipher_name} S-BOX")
    print("=" * 80)
    print(f"S-box: {sbox}\n")

    # --------------------------------------------------------------------------
    # Step 1: Compute DDT and Extract Valid Differential Transitions
    # --------------------------------------------------------------------------
    ddt, n = find_ddt(sbox)
    print_ddt(ddt, cipher_name)

    valid_points = extract_valid_transitions(ddt, n)
    total_possible = 2 ** (2 * n)
    print("=" * 80)
    print(f"STEP 1: VALID DDT TRANSITIONS ({len(valid_points)} / {total_possible} valid points)")
    print("=" * 80)
    for i, pt in enumerate(valid_points[:10], 1):
        in_part = tuple(pt[:n])
        out_part = tuple(pt[n:])
        print(f"Transition {i:3d}: {in_part} -> {out_part}  {pt}")
    if len(valid_points) > 10:
        print(f"... and {len(valid_points) - 10} more valid transitions.")
    print()

    # --------------------------------------------------------------------------
    # Step 2: Compute H-Representation (Convex Hull Inequalities via SciPy)
    # --------------------------------------------------------------------------
    raw_inequalities = compute_h_representation(valid_points)
    print("=" * 80)
    print(f"STEP 2: GENERATED H-REPRESENTATION INEQUALITIES ({len(raw_inequalities)} raw constraints)")
    print("=" * 80)
    
    # Differential variables naming convention: x for in_diff, y for out_diff
    var_names = [f"x{i}" for i in reversed(range(n))] + [f"y{i}" for i in reversed(range(n))]
    for i, ineq in enumerate(raw_inequalities[:5], 1):
        print(f"Facet {i:3d}: {format_algebraic(ineq, var_names=var_names)}")
    if len(raw_inequalities) > 5:
        print(f"... and {len(raw_inequalities) - 5} more facet inequalities.")
    print()

    # --------------------------------------------------------------------------
    # Step 3: Reduce Inequalities using Algorithm 1 (Greedy Elimination)
    # --------------------------------------------------------------------------
    reduced_ineqs = reduce_inequalities(valid_points, raw_inequalities)
    print("=" * 80)
    print(f"STEP 3: REDUCED LINEAR INEQUALITIES ({len(reduced_ineqs)} constraints)")
    print("=" * 80)
    for i, ineq in enumerate(reduced_ineqs, 1):
        print(f"Constraint {i:2d}: {format_algebraic(ineq, var_names=var_names)}")
    print("=" * 80)

    # --------------------------------------------------------------------------
    # Step 4: Verification of the Reduced Inequalities
    # --------------------------------------------------------------------------
    print("\nVERIFYING REDUCED INEQUALITIES...")
    passed, msg = verify_inequalities(valid_points, reduced_ineqs, n)
    if passed:
        print(f"SUCCESS: {msg}")
    else:
        print(f"FAILURE: {msg}")
    print("=" * 80)


if __name__ == "__main__":
    main()
