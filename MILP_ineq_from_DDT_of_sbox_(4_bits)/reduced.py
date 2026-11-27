"""
Inequality Size Reduction (Algorithm 1 from Sun et al., Paper 2016-857)
Reference: "Applying MILP to Division Property: An Approach for Division Trails"

This module implements Algorithm 1: greedy selection of a minimal subset of
linear inequalities from the H-representation such that its feasible binary
solutions restricted to {0, 1}^n are exactly the set of valid division trails.
"""

class Reduce:
    """
    Selects a minimal subset of linear inequalities from H-representation
    such that its feasible solutions in {0, 1}^n are exactly the division trails.
    """
    def __init__(self, divtrails, inequalities):
        self.divtrials = divtrails
        self.inequalities = inequalities

    @staticmethod
    def integer_to_bitlist(n, length):
        return [int(b) for b in format(n, f"0{length}b")]

    @staticmethod
    def value_of_expression(point, ineq):
        # ineq is (a0, a1, ..., a_{n-1}, c)
        # value is sum(a_i * point_i) + c
        return sum(x * y for x, y in zip(point, ineq[:-1])) + ineq[-1]

    def reduce(self):
        points = self.divtrials
        length = len(points[0])

        # apoints: all 2^length possible binary points
        apoints = [self.integer_to_bitlist(i, length) for i in range(2**length)]
        # cpoints: invalid division trails B = {0, 1}^n \ A
        points_set = set(tuple(p) for p in points)
        cpoints = [p for p in apoints if tuple(p) not in points_set]

        available_ineqs = [list(iq) for iq in self.inequalities]
        reduced_ineqs = []

        # Greedy selection: pick inequality that violates the maximum remaining invalid points
        while len(cpoints) > 0:
            best_violated_points = []
            best_ineq = None

            for l in available_ineqs:
                violated = [p for p in cpoints if self.value_of_expression(p, l) < 0]
                if len(violated) > len(best_violated_points):
                    best_violated_points = violated
                    best_ineq = l

            if best_ineq is None or len(best_violated_points) == 0:
                break

            violated_set = set(tuple(p) for p in best_violated_points)
            cpoints = [p for p in cpoints if tuple(p) not in violated_set]
            reduced_ineqs.append(tuple(best_ineq))
            available_ineqs.remove(best_ineq)

        return reduced_ineqs


def reduce_inequalities(vertices, inequalities):
    """
    Convenience function to run Algorithm 1 reduction on a set of inequalities.
    """
    reducer = Reduce(vertices, inequalities)
    return reducer.reduce()
