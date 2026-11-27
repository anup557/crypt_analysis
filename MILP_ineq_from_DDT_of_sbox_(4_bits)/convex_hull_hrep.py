"""
H-Representation of Convex Hull for Division Trails (Pure SciPy)
Reference: "Applying MILP to Division Property: An Approach for Division Trails" (Ling Sun et al., IACR ePrint 2016/857)

This module computes the H-representation (facet-defining linear inequalities)
from a set of vertices (division trails) using standard SciPy ConvexHull.
It does NOT require SageMath.
"""

import math
import numpy as np
from scipy.spatial import ConvexHull

def compute_h_representation(vertices):
    """
    Computes H-representation using standard SciPy ConvexHull.
    Normalizes each hyperplane equation into coprime integer coefficients:
        A * x + b >= 0
    Returns a sorted list of integer tuples: (a3, a2, a1, a0, b3, b2, b1, b0, b)
    """
    pts = np.array(vertices, dtype=float)
    hull = ConvexHull(pts)

    unique_ineqs = set()

    for eq in hull.equations:
        # SciPy's equation: normal . x + offset <= 0
        # In H-representation: A . x + b >= 0  =>  A = -normal, b = -offset
        r = -eq

        nz = [abs(x) for x in r if abs(x) > 1e-5]
        if not nz:
            continue

        m = min(nz)
        scaled = r / m

        # Find exact coprime integer representation
        for d in range(1, 100):
            test = scaled * d
            rounded = np.round(test)
            if np.max(np.abs(test - rounded)) < 1e-4:
                ints = rounded.astype(int)
                g = 0
                for v in ints:
                    g = math.gcd(g, abs(v))
                unique_ineqs.add(tuple(int(x) for x in (ints // g)))
                break

    return sorted(list(unique_ineqs))


def format_algebraic(ineq, var_names=None):
    """
    Format inequality as a standard MILP algebraic constraint:
        e.g., -4a3 - a2 - a1 - 4a0 + 2b3 + 2b2 + 2b1 - 3b0 + 7 >= 0
    """
    coeffs = ineq[:-1]
    const = ineq[-1]
    n = len(coeffs)
    if var_names is None:
        if n == 8:
            var_names = ['a3', 'a2', 'a1', 'a0', 'b3', 'b2', 'b1', 'b0']
        else:
            var_names = [f'x{i}' for i in range(n)]

    terms = []
    for c, v in zip(coeffs, var_names):
        if c == 0:
            continue
        elif c == 1:
            terms.append(f"+ {v}" if terms else f"{v}")
        elif c == -1:
            terms.append(f"- {v}" if terms else f"-{v}")
        elif c > 0:
            terms.append(f"+ {c}{v}" if terms else f"{c}{v}")
        else:
            terms.append(f"- {-c}{v}" if terms else f"-{-c}{v}")

    if const > 0:
        terms.append(f"+ {const}" if terms else f"{const}")
    elif const < 0:
        terms.append(f"- {-const}" if terms else f"-{-const}")

    lhs = " ".join(terms) if terms else "0"
    return f"{lhs} >= 0"
