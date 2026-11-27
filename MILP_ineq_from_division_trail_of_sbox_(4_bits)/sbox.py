"""
Division Trails Generation for S-boxes (Algorithm 2, Paper 2016-857)
Reference: "Applying MILP to Division Property: An Approach for Division Trails" (Ling Sun et al.)

This module implements Algorithm 2 to generate all valid division trails
from an S-box lookup table via ANF and bit-product properties.
"""

class Sbox:
    def __init__(self, sbox):
        self.sbox = sbox
        self.SBOXSIZE = self.SboxSize()

    def SboxSize(self):
        """
        Calculate the bit-size n of the S-box (length must be 2^n).
        """
        s = format(len(self.sbox), "b")
        num_of_1s = s.count("1")
        assert num_of_1s == 1, "S-box length must be a power of 2"
        return len(s) - 1

    def BitProduct(self, u, x):
        """
        Evaluates the bit-product function Pi_u(x).
        Returns 1 if (u & x) == u else 0.
        """
        return 1 if (u & x) == u else 0

    def GetTruthTable(self, u):
        """
        Retrieve truth table of boolean function Pi_u(S(x)).
        """
        temp = [u for _ in range(len(self.sbox))]
        table = list(map(self.BitProduct, temp, self.sbox))
        return table

    def ProcessTable(self, table):
        """
        Fast Moebius Transform to calculate the ANF from the truth table.
        """
        for i in range(0, self.SBOXSIZE):
            for j in range(0, 2**i):
                for k in range(0, 2**(self.SBOXSIZE - 1 - i)):
                    idx1 = k + 2**(self.SBOXSIZE - 1 - i) + j * (2**(self.SBOXSIZE - i))
                    idx2 = k + j * (2**(self.SBOXSIZE - i))
                    table[idx1] ^= table[idx2]

    def CreatANF(self):
        """
        Calculate the Algebraic Normal Form (ANF) of the coordinates
        and coordinate products of the S-box.
        """
        ANF = [[] for _ in range(len(self.sbox))]
        for i in range(1, len(self.sbox)):
            table = self.GetTruthTable(i)
            self.ProcessTable(table)
            sqr = [j for j in range(len(self.sbox)) if table[j] != 0]
            ANF[i] = sqr
        return ANF

    def CreateDivisionTrails(self):
        """
        Generate all valid division trails (Algorithm 2).
        Returns a list of 2n-bit binary lists [a_{n-1}, ..., a_0, b_{n-1}, ..., b_0].
        """
        ANF = self.CreatANF()
        INDP = []

        # Add zero vector [0, 0, ..., 0]
        sqr = [0 for _ in range(2 * self.SBOXSIZE)]
        INDP.append(sqr)

        # Start from non-zero vectors
        for i in range(1, len(self.sbox)):
            sqn = []
            for j in range(1, len(self.sbox)):
                flag = any((i | entry) == entry for entry in ANF[j])
                if flag:
                    sqn1 = []
                    flag_add = True
                    for t1 in sqn:
                        if (t1 | j) == j:
                            flag_add = False
                            break
                        elif (t1 | j) == t1:
                            sqn1.append(t1)
                    if flag_add:
                        for t2 in sqn1:
                            sqn.remove(t2)
                        sqn.append(j)

            for num in sqn:
                a = format(i, "0256b")
                b = format(num, "0256b")
                a = list(reversed(list(map(int, list(a)))))
                b = list(reversed(list(map(int, list(b)))))
                a = a[0:self.SBOXSIZE]
                b = b[0:self.SBOXSIZE]
                a.reverse()
                b.reverse()
                INDP.append((a + b))

        return INDP
