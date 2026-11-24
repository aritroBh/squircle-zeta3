lll_code = r'''
# lll_rel.py -- integer relation finding via flint fmpz_mat.lll (fast, exact)
import mpmath as mp
import flint

def find_relations(xs, prec_digits, max_report=6, scale_pow=None):
    n = len(xs)
    if scale_pow is None:
        scale_pow = prec_digits
    scale = mp.mpf(10) ** scale_pow
    B = []
    for i in range(n):
        row = [0] * (n + 1)
        row[i] = 1
        row[n] = int(mp.nint(xs[i] * scale))
        B.append(row)
    M = flint.fmpz_mat(B)
    R = M.lll()
    out = []
    rows = R.tolist()
    for row in rows[:max_report]:
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
