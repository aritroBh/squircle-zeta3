lll_code = r'''
# lll_rel.py -- integer relation finding via floating-point LLL (mpmath Gram-Schmidt,
# generous guard digits, re-orthogonalization)
import mpmath as mp

def lll_reduce_fp(Bint, prec, delta=mp.mpf('0.99')):
    with mp.workprec(prec):
        B = [list(map(int, b)) for b in Bint]
        n = len(B)
        dim = len(B[0])
        def gs():
            Bm = [[mp.mpf(x) for x in row] for row in B]
            Bstar = []; Bs = []; mu = [[mp.mpf(0)] * n for _ in range(n)]
            for i in range(n):
                bstar = list(Bm[i])
                for j in range(i):
                    if Bs[j]:
                        mu[i][j] = mp.fsum(Bm[i][t] * Bstar[j][t] for t in range(dim)) / Bs[j]
                    bstar = [x - mu[i][j] * y for x, y in zip(bstar, Bstar[j])]
                # re-orthogonalize once
                for j in range(i):
                    if Bs[j]:
                        corr = mp.fsum(bstar[t] * Bstar[j][t] for t in range(dim)) / Bs[j]
                    else:
                        corr = 0
                    if corr:
                        mu[i][j] += corr
                        bstar = [x - corr * y for x, y in zip(bstar, Bstar[j])]
                Bstar.append(bstar)
                Bs.append(mp.fsum(x * x for x in bstar))
            return mu, Bs
        k = 1
        while k < n:
            mu, Bs = gs()
            for j in range(k - 1, -1, -1):
                r = mp.nint(mu[k][j])
                if r:
                    ri = int(r)
                    B[k] = [x - ri * y for x, y in zip(B[k], B[j])]
                    mu, Bs = gs()
            if Bs[k] >= (delta - mu[k][k - 1] ** 2) * Bs[k - 1]:
                k += 1
            else:
                B[k], B[k - 1] = B[k - 1], B[k]
                k = max(k - 1, 1)
        return B

def find_relations(xs, prec_digits, max_report=6, guard=120):
    n = len(xs)
    scale = mp.mpf(10) ** prec_digits
    B = []
    for i in range(n):
        row = [0] * (n + 1)
        row[i] = 1
        row[n] = int(mp.nint(xs[i] * scale))
        B.append(row)
    R = lll_reduce_fp(B, prec_digits + guard)
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
l2 = mp.log(2); z3 = mp.zeta(3); p2 = mp.pi**2
xs = [mp.polylog(3, mp.mpf(1)/2), z3, p2*l2, l2**3, l2**2, mp.pi**4]
for a, resid in find_relations(xs, 90, 3):
    print(a, " resid:", mp.nstr(resid, 3))
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=300)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
