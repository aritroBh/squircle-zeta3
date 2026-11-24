
dixon_module = '''"""thx1_dixon.py -- exact rational solve of S h = b (S square invertible over Q, sparse)
via Dixon p-adic lifting with flint acceleration. Result verified exactly afterwards.
"""
from fractions import Fraction as Fr
import flint
from math import gcd
from functools import reduce

def lcm(a, b): return a // gcd(a, b) * b

def prep_system(Srows, b, r):
    """Srows: dict {row: {col: Fraction}}, b: list of Fraction (len r).
    Returns integer row-scaled system: A[i][j] ints, bb[i] ints, and scale."""
    A = {}
    bb = [0]*r
    for i in range(r):
        row = Srows.get(i, {})
        den = 1
        for v in row.values():
            den = lcm(den, v.denominator)
        den = lcm(den, b[i].denominator)
        A[i] = {j: int(v*den) for j, v in row.items()}
        bb[i] = int(b[i]*den)
    return A, bb

def dixon_solve(A, bb, r, p, nsteps_max=40, target_bits=None):
    """Solve A x = bb over Q. Returns (x_hat, p^k, x_list mod p^k) for reconstruction.
    A: {row: {col: int}}, bb: list of int."""
    # Ainv mod p via flint
    Af = flint.nmod_mat(r, r, p)
    for i, row in A.items():
        for j, v in row.items():
            Af[i, j] = v % p
    Ainv = Af.inv()
    # lifting
    x_hat = [0]*r       # current approximation (big ints)
    pk = 1
    bwork = list(bb)    # current residual/p^k
    xcols = []
    for step in range(nsteps_max):
        # x_i = Ainv * bwork mod p
        bv = flint.nmod_mat(r, 1, p)
        for i in range(r):
            bv[i, 0] = bwork[i] % p
        xi = (Ainv * bv)
        xil = [int(xi[i, 0]) % p for i in range(r)]
        # x_hat += xil * pk ; bwork = (bwork - A xil)/p
        for i in range(r):
            x_hat[i] += xil[i]*pk
        newb = [0]*r
        for i in range(r):
            s = bwork[i]
            for j, v in A[i].items():
                s -= v*xil[j]
            assert s % p == 0
            newb[i] = s//p
        bwork = newb
        pk *= p
        if target_bits and pk.bit_length() > 2*target_bits + 4:
            break
    return x_hat, pk
'''
open('thx1_dixon.py','w').write(dixon_module)
print("written thx1_dixon.py")
