t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
import sympy as sp, mpmath as mp, itertools
mp.mp.dps = 100

def log_matrix(c0, prec=90):
    xx = sp.Symbol('xx')
    mpoly = sp.minpoly(c0, xx)
    roots = sp.real_roots(sp.Poly(mpoly, xx))
    cvals = [mp.mpf(str(r.evalf(prec))) for r in roots]
    M = []
    for cv in cvals:
        M.append([mp.log(abs(v)) for v in (mp.mpf(2), cv, 1 - cv, 1 + cv)])
    return mpoly, M

for name, c0 in [("sin(pi/5)", sp.sin(sp.pi/5)), ("sin(pi/8)", sp.sin(sp.pi/8)),
                 ("sin(pi/12)", sp.sin(sp.pi/12)), ("cos(pi/12)", sp.cos(sp.pi/12)),
                 ("1/sqrt2", 1/sp.sqrt(2))]:
    mpoly, M = log_matrix(c0)
    d = len(M)
    print(f"=== {name}: deg {d}")
    # rank via all minors of max size
    maxminor = mp.mpf(0); bestk = 0
    for k in range(1, min(4, d) + 1):
        best = mp.mpf(0)
        for rows in itertools.combinations(range(d), k):
            for cols in itertools.combinations(range(4), k):
                sub = mp.matrix([[M[j][i] for i in cols] for j in rows])
                detv = mp.fabs(mp.det(sub))
                if detv > best: best = detv
        print(f"   best {k}x{k} minor |det| = {mp.nstr(best, 6)}")
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=300)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
