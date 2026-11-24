"""Finite differences in p (no differentiation under the integral): central 5th difference on
p = 1 + (k - 1/2) h, k = -2..3, Richardson in h^2.  Own scalar quadrature (direct5.scalar_P)."""
import sys, mpmath as mp, sympy as sp
import direct5 as D, formula5 as FM
from run_hex5 import CASES
name = sys.argv[1]
mp.mp.dps = 80
V, w = CASES[name]
H = D.setup([(FM.N(a), FM.N(b)) for a, b in V], [FM.N(x) for x in w])
def d5(h):
    pts = [mp.mpf(5) / 2, mp.mpf(3) / 2, mp.mpf(1) / 2, -mp.mpf(1) / 2, -mp.mpf(3) / 2, -mp.mpf(5) / 2]
    cf = [1, -5, 10, -10, 5, -1]
    return mp.fsum(c * D.scalar_P(H, 1 + s * h) for c, s in zip(cf, pts)) / h ** 5
hs = [mp.mpf(1) / 1000 / 2 ** k for k in range(6)]
T = [[d5(h)] for h in hs]
for i in range(1, len(hs)):
    for j in range(1, i + 1):
        T[i].append(T[i][j - 1] + (T[i][j - 1] - T[i - 1][j - 1]) / (4 ** j - 1))
ref = mp.mpf(sys.argv[2])
for i in range(len(hs)):
    print('h=%s  levels:' % mp.nstr(hs[i], 3), [mp.nstr(T[i][j] - ref, 3) for j in range(i + 1)])
print('%s: best Richardson value %s, minus jet value %s' % (name, mp.nstr(T[-1][-1], 30), mp.nstr(T[-1][-1] - ref, 3)))
