#!/usr/bin/env python3
"""
MILP Modeling for S-box Differential Distribution Tables (*-DDT) via Espresso
=============================================================================
Based on the paper:
  "MILP Modeling for (Large) S-boxes to Optimize Probability of Differential Characteristics"
  by Ahmed Abdelkhalek, Yu Sasaki, Yosuke Todo, Mohamed Tolba, Amr M. Youssef
  IACR Transactions on Symmetric Cryptology (ToSC), Vol. 2017, No. 4, pp. 99-129.

This module implements the end-to-end pipeline:
  1. Input: S-box lookup table S: F_2^n -> F_2^m
  2. Step 1: Compute the DDT (Differential Distribution Table)
  3. Step 2: Generate the *-DDT (Binary table: 1 if DDT[dx][dy] > 0, else 0)
  4. Step 3: Extract impossible differential minterms and construct the SOP
             (Sum-of-Products) of the complement function f_bar
  5. Step 4: Use Espresso logic minimization to find the minimal prime implicants
  6. Step 5: Apply De Morgan's duality to obtain the POS (Product-of-Sums / CNF) clauses of f
  7. Step 6: Convert POS clauses to MILP linear inequalities with coefficients in {-1, 0, 1}
  8. Extension (Section 4 of paper): Generate pb-DDT slices and Big-M conditional constraints
"""

import math
from typing import List, Tuple, Dict, Any, Optional

try:
    from pyeda.boolalg.espresso import espresso, FTYPE
    PYEDA_AVAILABLE = True
except ImportError:
    PYEDA_AVAILABLE = False


class SBoxMILPModeler:
    """
    Automated MILP constraint generator for cryptographic S-box Differential Distribution Tables.
    """

    def __init__(self, sbox: List[int], n_bits: Optional[int] = None, m_bits: Optional[int] = None):
        """
        :param sbox: Lookup table representing the S-box mapping x -> S[x]
        :param n_bits: Input bit width (inferred from len(sbox) if None)
        :param m_bits: Output bit width (inferred from max(sbox) if None)
        """
        self.sbox = list(sbox)
        self.n_in = n_bits if n_bits is not None else int(math.ceil(math.log2(len(sbox))))
        max_val = max(self.sbox) if self.sbox else 0
        self.n_out = m_bits if m_bits is not None else (int(math.ceil(math.log2(max_val + 1))) if max_val > 0 else self.n_in)

        # Variable names in standard cryptanalysis notation:
        # Input bits:  x_{n_in-1}, ..., x_0  (MSB to LSB)
        # Output bits: y_{n_out-1}, ..., y_0 (MSB to LSB)
        self.in_vars = [f"x{i}" for i in range(self.n_in - 1, -1, -1)]
        self.out_vars = [f"y{j}" for j in range(self.n_out - 1, -1, -1)]
        self.all_vars = self.in_vars + self.out_vars
        self.total_vars = self.n_in + self.n_out

        self.ddt: List[List[int]] = []
        self.star_ddt: List[List[int]] = []
        self.impossible_minterms: List[Tuple[int, int]] = []
        self.minimized_cubes: List[Tuple[Tuple[int, ...], Tuple[int]]] = []
        self.pos_clauses: List[str] = []
        self.inequalities: List[Dict[str, Any]] = []

    # -------------------------------------------------------------------------
    # Step 1: Differential Distribution Table (DDT)
    # -------------------------------------------------------------------------
    def compute_ddt(self) -> List[List[int]]:
        """
        Compute DDT[dx][dy] = # { x in F_2^n | S[x] ^ S[x ^ dx] == dy }
        """
        in_size = 1 << self.n_in
        out_size = 1 << self.n_out
        self.ddt = [[0] * out_size for _ in range(in_size)]

        for dx in range(in_size):
            for x in range(in_size):
                dy = self.sbox[x] ^ self.sbox[x ^ dx]
                self.ddt[dx][dy] += 1
        return self.ddt

    # -------------------------------------------------------------------------
    # Step 2: *-DDT (Binary Possibility Matrix)
    # -------------------------------------------------------------------------
    def compute_star_ddt(self) -> List[List[int]]:
        """
        Compute *-DDT[dx][dy] = 1 if DDT[dx][dy] > 0 else 0
        """
        if not self.ddt:
            self.compute_ddt()

        in_size = 1 << self.n_in
        out_size = 1 << self.n_out
        self.star_ddt = [
            [1 if self.ddt[dx][dy] > 0 else 0 for dy in range(out_size)]
            for dx in range(in_size)
        ]
        return self.star_ddt

    # -------------------------------------------------------------------------
    # Step 3: Extract Impossible Minterms & Formulate SOP
    # -------------------------------------------------------------------------
    def extract_impossible_minterms(self) -> List[Tuple[int, int]]:
        """
        Extract all (dx, dy) pairs where *-DDT[dx][dy] == 0.
        These are the minterms of f_bar (the complement function).
        """
        if not self.star_ddt:
            self.compute_star_ddt()

        in_size = 1 << self.n_in
        out_size = 1 << self.n_out
        self.impossible_minterms = [
            (dx, dy)
            for dx in range(in_size)
            for dy in range(out_size)
            if self.star_ddt[dx][dy] == 0
        ]
        return self.impossible_minterms

    def build_sop_terms(self) -> List[str]:
        """
        Generate symbolic SOP product terms for the impossible transitions:
          f_bar(x, y) = Or( And(l_1, l_2, ...), ... )
        """
        if not self.impossible_minterms:
            self.extract_impossible_minterms()

        sop_strings = []
        for dx, dy in self.impossible_minterms:
            literals = []
            for bit in range(self.n_in - 1, -1, -1):
                val = (dx >> bit) & 1
                var = f"x{bit}"
                literals.append(var if val == 1 else f"~{var}")
            for bit in range(self.n_out - 1, -1, -1):
                val = (dy >> bit) & 1
                var = f"y{bit}"
                literals.append(var if val == 1 else f"~{var}")
            sop_strings.append(" & ".join(literals))
        return sop_strings

    # -------------------------------------------------------------------------
    # Step 4: Logic Minimization via Espresso
    # -------------------------------------------------------------------------
    def minimize_via_espresso(self) -> List[Tuple[Tuple[int, ...], Tuple[int]]]:
        """
        Minimize the SOP of impossible transitions using the Espresso algorithm.
        Positional cube representation in PyEDA:
          1 = true (var == 1)
          2 = false (var == 0)
          3 = don't care (-)
        """
        if not PYEDA_AVAILABLE:
            raise RuntimeError("PyEDA library is required to run Espresso. Please install via: pip install pyeda")

        if not self.impossible_minterms:
            self.extract_impossible_minterms()

        cover = set()
        for dx, dy in self.impossible_minterms:
            invec = []
            for bit in range(self.n_in - 1, -1, -1):
                invec.append(1 if ((dx >> bit) & 1) else 2)
            for bit in range(self.n_out - 1, -1, -1):
                invec.append(1 if ((dy >> bit) & 1) else 2)
            cover.add((tuple(invec), (1,)))

        self.minimized_cubes = list(espresso(self.total_vars, 1, cover, intype=FTYPE))
        return self.minimized_cubes

    # -------------------------------------------------------------------------
    # Step 5 & 6: De Morgan POS Clauses -> Linear Inequalities
    # -------------------------------------------------------------------------
    def generate_pos_and_inequalities(self) -> Tuple[List[str], List[Dict[str, Any]]]:
        """
        Converts minimized SOP implicants into POS clauses and linear inequalities:
          1. Negate each product term using De Morgan's Law:
             ~ (And_{i in I_1} v_i & And_{j in I_0} ~v_j) = Or_{i in I_1} ~v_i | Or_{j in I_0} v_j
          2. Convert to linear inequality:
             sum_{j in I_0} v_j + sum_{i in I_1} (1 - v_i) >= 1
             ==> sum_{j in I_0} v_j - sum_{i in I_1} v_i >= 1 - |I_1|
        """
        if not self.minimized_cubes:
            self.minimize_via_espresso()

        self.pos_clauses = []
        self.inequalities = []

        for in_part, out_part in self.minimized_cubes:
            clause_lits = []
            pos_vars = []  # +1 coefficient (var was 0 in impossible cube)
            neg_vars = []  # -1 coefficient (var was 1 in impossible cube)

            for idx, val in enumerate(in_part):
                var = self.all_vars[idx]
                if val == 1:
                    # var was 1 in impossible cube -> ~var in POS -> (1 - var) in inequality
                    clause_lits.append(f"~{var}")
                    neg_vars.append(var)
                elif val == 2:
                    # var was 0 in impossible cube -> var in POS -> var in inequality
                    clause_lits.append(var)
                    pos_vars.append(var)
                # val == 3 is don't care (-), omitted

            rhs = 1 - len(neg_vars)

            # Format linear string
            terms = [f"+ {v}" for v in pos_vars] + [f"- {v}" for v in neg_vars]
            expr_str = " ".join(terms)
            if expr_str.startswith("+ "):
                expr_str = expr_str[2:]
            elif expr_str.startswith("- "):
                expr_str = "-" + expr_str[2:]

            ineq_str = f"{expr_str} >= {rhs}"
            pos_str = "(" + " | ".join(clause_lits) + ")"

            self.pos_clauses.append(pos_str)
            self.inequalities.append({
                "inequality_str": ineq_str,
                "pos_vars": pos_vars,
                "neg_vars": neg_vars,
                "rhs": rhs,
                "pos_clause": pos_str
            })

        return self.pos_clauses, self.inequalities

    # -------------------------------------------------------------------------
    # Extension (Section 4): pb-DDT Decomposition with Big-M Constraints
    # -------------------------------------------------------------------------
    def generate_pb_inequalities(self, M: Optional[int] = None) -> Dict[int, Dict[str, Any]]:
        """
        Decomposes DDT into probability-specific slices (pb-DDTs).
        Generates conditional Big-M inequalities:
          <a, (x, y)> + M * (1 - Q_pb) >= rhs
        where Q_pb in {0, 1} is an indicator variable for that probability.
        """
        if not self.ddt:
            self.compute_ddt()

        if M is None:
            M = self.total_vars  # M = 2n is sufficient and tight

        # Group distinct non-zero transition counts
        counts = set()
        for row in self.ddt:
            for val in row:
                if val > 0:
                    counts.add(val)

        in_size = 1 << self.n_in
        out_size = 1 << self.n_out
        total_pairs = in_size  # S-box input size

        pb_results = {}

        for count in sorted(counts):
            prob_float = count / total_pairs
            log2_p = math.log2(prob_float)
            q_var = f"Q_p_{count}"

            # Construct pb-DDT where entry is 1 iff ddt[dx][dy] == count
            cover = set()
            for dx in range(in_size):
                for dy in range(out_size):
                    if self.ddt[dx][dy] != count:
                        invec = []
                        for bit in range(self.n_in - 1, -1, -1):
                            invec.append(1 if ((dx >> bit) & 1) else 2)
                        for bit in range(self.n_out - 1, -1, -1):
                            invec.append(1 if ((dy >> bit) & 1) else 2)
                        cover.add((tuple(invec), (1,)))

            cubes = list(espresso(self.total_vars, 1, cover, intype=FTYPE))
            slice_inequalities = []

            for in_part, _ in cubes:
                pos_vars = []
                neg_vars = []
                for idx, val in enumerate(in_part):
                    var = self.all_vars[idx]
                    if val == 1:
                        neg_vars.append(var)
                    elif val == 2:
                        pos_vars.append(var)

                rhs = 1 - len(neg_vars)
                # Conditional inequality: sum(pos) - sum(neg) - M * Q_pb >= rhs - M
                terms = [f"+ {v}" for v in pos_vars] + [f"- {v}" for v in neg_vars] + [f"+ {M} {q_var}"]
                expr_str = " ".join(terms)
                if expr_str.startswith("+ "):
                    expr_str = expr_str[2:]
                ineq_str = f"{expr_str} >= {rhs + M}"

                slice_inequalities.append({
                    "inequality_str": ineq_str,
                    "pos_vars": pos_vars,
                    "neg_vars": neg_vars,
                    "q_var": q_var,
                    "M": M,
                    "rhs": rhs
                })

            pb_results[count] = {
                "count": count,
                "probability": prob_float,
                "log2_p": log2_p,
                "q_var": q_var,
                "num_constraints": len(slice_inequalities),
                "constraints": slice_inequalities
            }

        return pb_results

    # -------------------------------------------------------------------------
    # Verification
    # -------------------------------------------------------------------------
    def verify_inequalities(self) -> Dict[str, Any]:
        """
        Mathematically verifies soundness and completeness:
          - 100% of possible transitions (*-DDT == 1) satisfy all inequalities.
          - 100% of impossible transitions (*-DDT == 0) violate at least one inequality.
        """
        if not self.inequalities:
            self.generate_pos_and_inequalities()

        in_size = 1 << self.n_in
        out_size = 1 << self.n_out

        valid_passes = 0
        valid_total = 0
        invalid_rejected = 0
        invalid_total = 0

        for dx in range(in_size):
            dx_bits = [(dx >> bit) & 1 for bit in range(self.n_in - 1, -1, -1)]
            for dy in range(out_size):
                dy_bits = [(dy >> bit) & 1 for bit in range(self.n_out - 1, -1, -1)]
                val_dict = {}
                for v, b in zip(self.in_vars, dx_bits):
                    val_dict[v] = b
                for v, b in zip(self.out_vars, dy_bits):
                    val_dict[v] = b

                all_sat = True
                for ineq in self.inequalities:
                    lhs = sum(val_dict[v] for v in ineq["pos_vars"]) - sum(val_dict[v] for v in ineq["neg_vars"])
                    if lhs < ineq["rhs"]:
                        all_sat = False
                        break

                if self.star_ddt[dx][dy] == 1:
                    valid_total += 1
                    if all_sat:
                        valid_passes += 1
                else:
                    invalid_total += 1
                    if not all_sat:
                        invalid_rejected += 1

        return {
            "valid_total": valid_total,
            "valid_satisfied": valid_passes,
            "all_valid_pass": (valid_total == valid_passes),
            "invalid_total": invalid_total,
            "invalid_rejected": invalid_rejected,
            "all_invalid_rejected": (invalid_total == invalid_rejected)
        }


# =============================================================================
# Demo & Execution
# =============================================================================
if __name__ == "__main__":
    import sys

    print("=" * 80)
    print("DEMO 1: 3-bit Illustrative S-box from Table 1 of ToSC 2017 Paper")
    print("=" * 80)
    sbox_3bit = [0x5, 0x3, 0x4, 0x6, 0x2, 0x7, 0x0, 0x1]
    modeler_3bit = SBoxMILPModeler(sbox_3bit)

    # 1. DDT
    modeler_3bit.compute_ddt()
    print(f"[*] DDT computed (8 x 8 table)")

    # 2. *-DDT
    star_ddt = modeler_3bit.compute_star_ddt()
    num_possible = sum(sum(r) for r in star_ddt)
    num_impossible = 64 - num_possible
    print(f"[*] *-DDT computed: {num_possible} possible, {num_impossible} impossible transitions")

    # 3. SOP terms
    sop = modeler_3bit.build_sop_terms()
    print(f"[*] Constructed SOP from {len(sop)} impossible minterms")
    print(f"    Example minterm: {sop[0]}")

    # 4 & 5 & 6. Espresso -> POS -> Linear Inequalities
    pos_clauses, inequalities = modeler_3bit.generate_pos_and_inequalities()
    print(f"[*] Espresso optimized {len(sop)} minterms -> {len(inequalities)} prime implicants / inequalities")
    print("\nGenerated POS Clauses and MILP Inequalities:")
    for i, ineq in enumerate(inequalities):
        print(f"  [{i+1:2d}] POS: {ineq['pos_clause']}")
        print(f"       MILP: {ineq['inequality_str']}")

    # Verification
    v_res = modeler_3bit.verify_inequalities()
    print(f"\nVerification against complete truth table:")
    print(f"  - Valid transitions satisfied: {v_res['valid_satisfied']}/{v_res['valid_total']} ({v_res['all_valid_pass']})")
    print(f"  - Impossible transitions rejected: {v_res['invalid_rejected']}/{v_res['invalid_total']} ({v_res['all_invalid_rejected']})")

    # Section 4 Demo: pb-DDT Slices
    print("\n" + "=" * 80)
    print("DEMO 2: Probability Slices (pb-DDT) with Big-M Conditional Constraints")
    print("=" * 80)
    pb_slices = modeler_3bit.generate_pb_inequalities()
    for count, data in pb_slices.items():
        if count == 8:
            continue  # ignore zero-difference transition (handled via Q)
        print(f"  - Slice count={count} (p = 2^{data['log2_p']:.1f}): {data['num_constraints']} inequalities, variable: {data['q_var']}")
    print("=" * 80)
