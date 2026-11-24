lll_code = r'''
# lll_rel.py -- exact integer LLL + integer relation finding
from fractions import Fraction
import math

def lll_reduce(B, delta=Fraction(3, 4)):
    """B: list of rows (lists of ints/Fractions). Returns LLL-reduced basis."""
    B = [list(map(Fraction, b)) for b in B]
    n = len(B)
    def dot(u, v): return sum(a * b for a, b in zip(u, v))
    def gs():
        Bs = []
        mu = [[Fraction(0)] * n for _ in range(n)]
        Bstar = []
        for i in range(n):
            bstar = list(B[i])
            for j in range(i):
                mu[i][j] = dot(B[i], Bstar[j]) / Bs[j]
                bstar = [x - mu[i][j] * y for x, y in zip(bstar, Bstar[j])]
            Bstar.append(bstar)
            Bs.append(dot(bstar, bstar))
        return mu, Bstar, Bs
    k = 1
    while k < n:
        mu, Bstar, Bs = gs()
        for j in range(k - 1, -1, -1):
            q = mu[k][j]
            r = q.numerator // q.denominator
            if 2 * abs(q - r) > 1:
                r += 1 if q > 0 else -1
            if q - r > Fraction(1, 2):
                r += 1
            if q - r < Fraction(-1, 2):
                r -= 1
            if r != 0:
                B[k] = [x - r * y for x, y in zip(B[k], B[j])]
                mu, Bstar, Bs = gs()
        if Bs[k] >= (delta - mu[k][k - 1] ** 2) * Bs[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            k = max(k - 1, 1)
    return B

def find_relations(xs, prec_digits, max_report=4):
    """xs: list of mpf. Returns short near-relations: list of (int vector, residual)."""
    import mpmath as mp
    n = len(xs)
    scale = mp.mpf(10) ** prec_digits
    B = []
    for i in range(n):
        row = [0] * (n + 1)
        row[i] = 1
        row[n] = int(mp.nint(xs[i] * scale))
        B.append(row)
    R = lll_reduce(B)
    out = []
    for row in R[:max_report]:
        a = row[:n]
        resid = mp.fsum(mp.mpf(int(c)) * x for c, x in zip(a, xs))
        out.append(([int(c) for c in a], resid, row[n]))
    return out
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/lll_rel.py","w").write(lll_code)

t = r'''
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
from lll_rel import find_relations
mp.mp.dps = 200
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O2 = Omega(mp.mpf(1)/2)
C2 = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
D2 = O2 - C2
print("D2 =", mp.nstr(D2, 40))
l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3); p2 = mp.pi**2
L2b = mp.polylog(2, mp.mpf(1)/3)
L3half = mp.polylog(3, mp.mpf(1)/2); L3b = mp.polylog(3, mp.mpf(1)/3)
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b,
         L3half*l2, L3b*l2, L3half*l3, L3b*l3]
print("LLL relations for [D2, products] at 180 digits:")
for a, resid, last in find_relations([D2] + basis, 180, 8):
    print("   vec[:4]... resid:", mp.nstr(resid, 3), " a0 =", a[0], " small coeffs:", max(abs(c) for c in a) < 1000)
    if resid == 0 or abs(resid) < mp.mpf('1e-150'):
        print("     RELATION:", a)
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_05e.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05e.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
