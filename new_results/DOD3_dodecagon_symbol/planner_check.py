# Planner check of DOD3: L(3,chi_12), L(2,chi_12) by an independent route (Hurwitz zeta), and the S-unit solution set over Q(sqrt3).
from mpmath import mp, mpf, zeta, nstr
import pickle, sympy as sp
mp.dps = 40
chi = {1: 1, 5: -1, 7: -1, 11: 1}
L = lambda s: sum(c*zeta(s, mpf(a)/12) for a, c in chi.items())/mpf(12)**s
print("L(3,chi12) =", nstr(L(3), 25), " L(2,chi12) =", nstr(L(2), 25))
sols = pickle.load(open('solsF_B8.pkl', 'rb'))
print("type", type(sols), "len", len(sols))
s3 = sp.sqrt(3)
ex = list(sols)[:3] if not isinstance(sols, dict) else list(sols.items())[:3]
print("sample", ex)
from fractions import Fraction as Fr
def smooth(q):
    q = abs(q)
    for n in (q.numerator, q.denominator):
        for p in (2, 3):
            while n % p == 0: n //= p
        if n != 1: return False
    return True
keys = set()
bad = 0
for (a, b, c) in sols:
    z = (Fr(a, c), Fr(b, c))                      # z = x + y*sqrt3
    w = (1 - z[0], -z[1])
    nz = z[0]**2 - 3*z[1]**2; nw = w[0]**2 - 3*w[1]**2
    if nz == 0 or nw == 0 or not (smooth(nz) and smooth(nw)): bad += 1
    keys.add(z)
print("all 129: z and 1-z are S-units (norms 2,3-smooth):", bad == 0)
print("closed under z -> 1-z:", all(((1 - x), -y) in keys for x, y in keys))
def inv(x, y):
    n = x*x - 3*y*y; return (x/n, -y/n)
print("closed under z -> 1/z:", all(inv(x, y) in keys for x, y in keys))
