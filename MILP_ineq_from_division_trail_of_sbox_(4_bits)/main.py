"""
End-to-End Division Trail Modeling Pipeline (Paper 2016-857)
Reference: "Applying MILP to Division Property: An Approach for Division Trails" (Ling Sun et al.)

Pipeline:
    1. Input S-box lookup table.
    2. Generate and print all Division Trails of the S-box (Algorithm 2).
    3. Generate and print all H-Representation linear inequalities (Algebraic MILP format).
    4. Generate and print the Reduced Inequalities (Algebraic MILP format via Algorithm 1).
"""

from sbox import Sbox
from convex_hull_hrep import compute_h_representation, format_algebraic
from reduced import reduce_inequalities

def main():
    # --------------------------------------------------------------------------
    # Input: S-box lookup table
    # --------------------------------------------------------------------------
    # cipher_name = "PRESENT"
    # sbox = [0xc, 0x5, 0x6, 0xb, 0x9, 0x0, 0xa, 0xd, 0x3, 0xe, 0xf, 0x8, 0x4, 0x7, 0x1, 0x2]

    # cipher_name = "GIFT"
    # sbox = [0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe]

    cipher_name = "BAKSHEESH"
    sbox = [0x3, 0x0, 0x6, 0xD, 0xB, 0x5, 0x8, 0xE, 0xC, 0xF, 0x9, 0x2, 0x4, 0xA, 0x7, 0x1]

    print("=" * 80)
    print(f"DIVISION TRAIL & INEQUALITY MODELING PIPELINE FOR {cipher_name} S-BOX")
    print("=" * 80)
    print(f"S-box lookup table: {sbox}\n")

    # --------------------------------------------------------------------------
    # Step 1: Generate Division Trails from S-box (Algorithm 2)
    # --------------------------------------------------------------------------
    sb = Sbox(sbox)
    sbox_size = sb.SBOXSIZE
    division_trails = sb.CreateDivisionTrails()

    print("=" * 80)
    print(f"STEP 1: DIVISION TRAILS OF THE S-BOX ({len(division_trails)} trails)")
    print("=" * 80)
    for i, trail in enumerate(division_trails, 1):
        a_part = tuple(trail[:sbox_size])
        b_part = tuple(trail[sbox_size:])
        print(f"Trail {i:3d}: {a_part} -> {b_part}  {trail}")
    print()

    # --------------------------------------------------------------------------
    # Step 2: Generate H-Representation (Linear Inequalities) of the Convex Hull
    # --------------------------------------------------------------------------
    inequalities = compute_h_representation(division_trails)

    print("=" * 80)
    print(f"STEP 2: GENERATED H-REPRESENTATION INEQUALITIES ({len(inequalities)} constraints)")
    print("=" * 80)
    for i, ineq in enumerate(inequalities, 1):
        print(f"Constraint {i:3d}: {format_algebraic(ineq)}")
    print()

    # --------------------------------------------------------------------------
    # Step 3: Reduce Inequalities using Algorithm 1 (Greedy Elimination)
    # --------------------------------------------------------------------------
    reduced_ineqs = reduce_inequalities(division_trails, inequalities)

    print("=" * 80)
    print(f"STEP 3: REDUCED INEQUALITIES ({len(reduced_ineqs)} constraints)")
    print("=" * 80)
    for i, ineq in enumerate(reduced_ineqs, 1):
        print(f"Constraint {i:2d}: {format_algebraic(ineq)}")
    print("=" * 80)

if __name__ == "__main__":
    main()
