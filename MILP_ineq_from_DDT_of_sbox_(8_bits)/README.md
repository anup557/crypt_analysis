# MILP Modeling for (Large) S-Boxes to Optimize Differential Characteristics

An automated Python framework to model the Differential Distribution Tables (DDT / \*-DDT) of arbitrary $n$-bit S-boxes as linear constraints for Mixed Integer Linear Programming (MILP) solvers.

This implementation is based on the research paper:
> **MILP Modeling for (Large) S-boxes to Optimize Probability of Differential Characteristics**  
> *Ahmed Abdelkhalek, Yu Sasaki, Yosuke Todo, Mohamed Tolba, and Amr M. Youssef*  
> IACR Transactions on Symmetric Cryptology (ToSC), Vol. 2017, No. 4, pp. 99–129  
> DOI: [10.13154/tosc.v2017.i4.99-129](https://doi.org/10.13154/tosc.v2017.i4.99-129)  
> Paper PDF: [`refs/lstennes,+ToSC2017_4_05.pdf`](refs/lstennes,+ToSC2017_4_05.pdf)

---

## Table of Contents
- [Overview](#overview)
- [How It Works](#how-it-works)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
  - [1. Running the Built-in Demo](#1-running-the-built-in-demo)
  - [2. Python API Usage](#2-python-api-usage)
  - [3. Probability Slicing ($p_b$-DDT) with Big-M Constraints](#3-probability-slicing-p_b-ddt-with-big-m-constraints)
  - [4. Verification of Inequalities](#4-verification-of-inequalities)
- [Precomputed AES Inequalities](#precomputed-aes-inequalities)
- [Empirical Benchmarks](#empirical-benchmarks)
- [Citation](#citation)

---

## Overview

In automated differential cryptanalysis using MILP:
1. **Word-wise (truncated) models** fail to capture specific S-box differential behavior, yielding false-positive trails that cannot be instantiated with actual differences.
2. **Convex hull H-representation** requires facet enumeration in $\mathbb{R}^{2n}$, which becomes computationally infeasible when S-box dimension $n \ge 5$ (e.g., 8-bit S-boxes like AES or SKINNY-128).
3. **Naive logical condition modeling** produces one inequality per impossible transition (up to $2^{16} = 65,536$ constraints for an 8-bit S-box), overloading solvers.

**The Solution:**  
This project maps the logical condition modeling problem to **Boolean logic minimization**. By formulating the impossible differential transitions as an unoptimized canonical Product-of-Sums (POS / CNF) and applying the **Espresso algorithm**, the number of constraints is reduced by orders of magnitude while strictly preserving 100% soundness and completeness.

---

## How It Works

```
   ┌─────────────────────────────────────────────────────────────────┐
   │ 1. S-box Lookup Table S: F_2^n -> F_2^m                         │
   └────────────────────────────────┬────────────────────────────────┘
                                    │
                                    ▼
   ┌─────────────────────────────────────────────────────────────────┐
   │ 2. Compute DDT & *-DDT (Binary table: 1 = Possible, 0 = Imposs) │
   └────────────────────────────────┬────────────────────────────────┘
                                    │
                                    ▼
   ┌─────────────────────────────────────────────────────────────────┐
   │ 3. Extract Impossible Minterms: Formulate SOP of complement f̄   │
   └────────────────────────────────┬────────────────────────────────┘
                                    │
                                    ▼
   ┌─────────────────────────────────────────────────────────────────┐
   │ 4. Logic Minimization via Espresso:                             │
   │    Reduces impossible minterms to minimal prime implicants P_k  │
   └────────────────────────────────┬────────────────────────────────┘
                                    │
                                    ▼
   ┌─────────────────────────────────────────────────────────────────┐
   │ 5. De Morgan's Duality: Negate P_k to obtain POS clauses of f   │
   │    f = (~P_1) ∧ (~P_2) ∧ ... ∧ (~P_K)                           │
   └────────────────────────────────┬────────────────────────────────┘
                                    │
                                    ▼
   ┌─────────────────────────────────────────────────────────────────┐
   │ 6. Linear Inequality Mapping (Coefficients in {-1, 0, 1}):      │
   │    ∑_{v ∈ Pos} v - ∑_{v ∈ Neg} v >= 1 - |Neg|                   │
   └─────────────────────────────────────────────────────────────────┘
```

### Why Minimize the Impossible Transitions?
- $f(\vec{x}, \vec{y}) = 1$ for **possible** transitions, and $0$ for **impossible** transitions.
- An MILP solver must enforce $f(\vec{x}, \vec{y}) = 1$.
- An AND-of-ORs (POS / CNF): $C_1 \land C_2 \land \dots \land C_K = 1 \iff C_1 \ge 1, C_2 \ge 1, \dots, C_K \ge 1$. Every clause becomes an independent linear inequality **without auxiliary binary variables**.
- Espresso is an SOP minimizer. Passing the **impossible transitions** (where $\overline{f} = 1$) to Espresso yields minimal product terms $P_k$. Negating each term via De Morgan's Law ($\overline{P_k}$) yields the exact minimal POS clauses of valid transitions!

---

## Repository Structure

| File | Description |
| :--- | :--- |
| [`sbox_to_milp.py`](sbox_to_milp.py) | Main Python module implementing DDT generation, logic minimization, POS conversion, and Big-M constraints. |
| [`aes_sbox_inequalities.txt`](aes_sbox_inequalities.txt) | Complete set of **8,269 precomputed linear inequalities** for the 8-bit AES S-box. |
| [`refs/lstennes,+ToSC2017_4_05.pdf`](refs/lstennes,+ToSC2017_4_05.pdf) | Original ToSC 2017 paper by Abdelkhalek et al. |
| `README.md` | Documentation and usage instructions. |

---

## Installation

The tool requires Python 3.8+ and the `pyeda` library (which includes the Berkeley Espresso C engine).

```bash
# Standard installation
pip install pyeda

# If building on modern GCC (GCC 14/15+)
CFLAGS="-Wno-incompatible-pointer-types" pip install pyeda
```

---

## Usage Guide

### 1. Running the Built-in Demo
Execute the script directly to run test benchmarks on the paper's 3-bit illustrative S-box:

```bash
python3 sbox_to_milp.py
```

Output:
```text
================================================================================
DEMO 1: 3-bit Illustrative S-box from Table 1 of ToSC 2017 Paper
================================================================================
[*] DDT computed (8 x 8 table)
[*] *-DDT computed: 23 possible, 41 impossible transitions
[*] Constructed SOP from 41 impossible minterms
[*] Espresso optimized 41 minterms -> 13 prime implicants / inequalities

Generated POS Clauses and MILP Inequalities:
  [ 1] POS: (x2 | x1 | ~y1 | ~y0)
       MILP: x2 + x1 - y1 - y0 >= -1
  [ 2] POS: (x1 | x0 | ~y1 | y0)
       MILP: x1 + x0 + y0 - y1 >= 0
  ...
Verification against complete truth table:
  - Valid transitions satisfied: 23/23 (True)
  - Impossible transitions rejected: 41/41 (True)
```

---

### 2. Python API Usage

The class `SBoxMILPModeler` can process any square ($n \times n$) or non-square ($n \times m$) S-box:

```python
from sbox_to_milp import SBoxMILPModeler

# Example: 4-bit PRESENT S-box
sbox_present = [0xc, 0x5, 0x6, 0xb, 0x9, 0x0, 0xa, 0xd, 0x3, 0xe, 0xf, 0x8, 0x4, 0x7, 0x1, 0x2]
modeler = SBoxMILPModeler(sbox_present, n_bits=4, m_bits=4)

# Generate POS clauses and linear inequalities
pos_clauses, inequalities = modeler.generate_pos_and_inequalities()

print(f"Generated {len(inequalities)} linear inequalities.")
for ineq in inequalities[:5]:
    print(ineq["inequality_str"])
```

#### Non-Square S-boxes ($n \ne m$)
Non-square S-boxes (such as the $6 \times 4$ DES S-boxes) are fully supported:

```python
# DES S-box S1 (6-bit input, 4-bit output)
des_s1 = [
    14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7,
    0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8,
    4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0,
    15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13
]
modeler = SBoxMILPModeler(des_s1, n_bits=6, m_bits=4)
pos_clauses, inequalities = modeler.generate_pos_and_inequalities()
print(f"DES S1: {len(inequalities)} inequalities generated.")
```

---

### 3. Probability Slicing ($p_b$-DDT) with Big-M Constraints

To optimize exact differential characteristic probabilities (Section 4 of the paper), the DDT is decomposed into probability slices. Each slice is guarded by an indicator variable $Q_{p_b} \in \{0, 1\}$ using Big-M conditional constraints:

$$\langle \vec{a}, (\vec{x}, \vec{y}) \rangle + M \cdot Q_{p_b} \ge \text{RHS} + M$$

```python
# Generate probability-sliced inequalities
pb_data = modeler.generate_pb_inequalities(M=16)

for count, data in pb_data.items():
    print(f"Count: {count} | Prob: 2^{data['log2_p']:.2f} | Constraints: {data['num_constraints']} | Indicator: {data['q_var']}")
```

---

### 4. Verification of Inequalities

The modeler includes a built-in mathematical verification function to guarantee 100% soundness (no false rejections of valid patterns) and completeness (no false acceptances of invalid patterns):

```python
results = modeler.verify_inequalities()
print("All valid patterns satisfied:", results["all_valid_pass"])
print("All impossible patterns rejected:", results["all_invalid_rejected"])
```

---

## Precomputed AES Inequalities

The file [`aes_sbox_inequalities.txt`](aes_sbox_inequalities.txt) contains the complete system of **8,269 inequalities** describing the differential propagation of the 8-bit AES S-box.

- **Variables:** Input bits $x_7, \dots, x_0$, output bits $y_7, \dots, y_0$.
- **Format:** Compatible with standard MILP solvers (Gurobi, SCIP, CPLEX, COIN-OR CBC).

Example:
```text
x5 + x1 + y6 + y2 - x7 - x6 - x3 - x2 - x0 - y7 - y4 - y3 - y0 >= -8
x6 + x4 + x3 + y6 + y5 + y2 + y1 + y0 - x5 - x2 - x1 - x0 - y7 - y3 >= -5
x7 + x2 + x0 + y6 + y3 + y1 + y0 - x6 - x5 - x4 - x1 - y7 - y4 - y2 >= -6
```

---

## Empirical Benchmarks

| S-box | Dimensions | Impossible Entries | Naive Inequalities | Espresso Inequalities | Soundness & Completeness |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Illustrative Toy** | $3 \times 3$ | 41 | 41 | **13** *(exact match to paper)* | 100% Verified |
| **PRESENT** | $4 \times 4$ | 159 | 159 | **37** | 100% Verified |
| **ASCON** | $5 \times 5$ | 694 | 694 | **60** | 100% Verified |
| **DES S1** | $6 \times 4$ | 448 | 448 | **122** | 100% Verified |
| **SKINNY-128** | $8 \times 8$ | 54,067 | 54,067 | **271** | 100% Verified |
| **AES** | $8 \times 8$ | 33,150 | 33,150 | **8,269** *(paper: 8,302)* | 100% Verified |

---

## Citation

If you use this code in academic research, please cite the original paper:

```bibtex
@article{ToSC:ASTTY17,
  author    = {Ahmed Abdelkhalek and
               Yu Sasaki and
               Yosuke Todo and
               Mohamed Tolba and
               Amr M. Youssef},
  title     = {{MILP} Modeling for (Large) S-boxes to Optimize Probability of Differential
               Characteristics},
  journal   = {{IACR} Trans. Symmetric Cryptol.},
  volume    = {2017},
  number    = {4},
  pages     = {99--129},
  year      = {2017},
  doi       = {10.13154/tosc.v2017.i4.99-129}
}
```
