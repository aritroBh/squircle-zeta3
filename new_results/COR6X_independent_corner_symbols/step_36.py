lll_code = r'''
# lll_rel.py -- exact integer LLL + integer relation finding
from fractions import Fraction

def lll_reduce(B, delta=Fraction(3, 4)):
    B = [list(map(Fraction, b)) for b in B]
    n = len(B)
    def dot(u, v): return sum(a * b for a, b in zip(u, v))
    def gs():
        Bstar = []; Bs = []; mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            bstar = list(B[i])
            for j in range(i):
                mu[i][j] = dot(B[i], Bstar[j]) / Bs[j] if Bs[j] else Fraction(0)
                bstar = [x - mu[i][j] * y for x, y in zip(bstar, Bstar[j])]
            Bstar.append(bstar); Bs.append(dot(bstar, bstar))
        return mu, Bstar, Bs
    k = 1
    while k < n:
        mu, Bstar, Bs = gs()
        for j in range(k - 1, -1, -1):
            if not Bs[j]:
                continue
            r = mu[k][j]
            fl = r.numerator // r.denominator
            best = fl + 1 if 2 * (r - fl) > 1 else fl
            if best:
                B[k] = [x - best * y for x, y in zip(B[k], B[j])]
                mu, Bstar, Bs = gs()
        if Bs[k] >= (delta - mu[k][k - 1] ** 2) * Bs[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            k = max(k - 1, 1)
    return B

def find_relations(xs, prec_digits, max_report=6):
    """short near-relations among reals xs (mpf list)."""
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
        a = [int(c) for c in row[:n]]
        resid = mp.fsum(mp.mpf(int(c)) * x for c, x in zip(a, xs))
        out.append((a, resid))
    return out
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/lll_rel.py","w").write(lll_code)

t = r'''
import sys; sys.path.insert(0, '.')
import mpmath as mp
from lll_rel import find_relations
mp.mp.dps = 100
# validation: Li_3(1/2) = 7/8 zeta(3) - pi^2/12 log2 + log^3 2 / 6
l2 = mp.log(2); z3 = mp.zeta(3); p2 = mp.pi**2
xs = [mp.polylog(3, mp.mpf(1)/2), z3, p2*l2, l2**3, l2**2, mp.pi**4]
for a, resid in find_relations(xs, 90, 3):
    print(a, " resid:", mp.nstr(resid, 3))
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
