t = r'''
import sympy as sp, mpmath as mp
mp.mp.dps = 120
c0 = sp.sin(sp.pi/5)
xx = sp.Symbol('xx')
mpoly = sp.minpoly(c0, xx)
roots = sp.real_roots(sp.Poly(mpoly, xx))
cvals = [mp.mpf(str(r.evalf(110))) for r in roots]
M = mp.matrix([[mp.log(abs(v)) for v in (mp.mpf(2), cv, 1-cv, 1+cv)] for cv in cvals])
print("roots:", [mp.nstr(c, 8) for c in cvals])
det = mp.det(M)
print("det =", mp.nstr(det, 10))
# kernel via svd-ish: null space
import numpy as np
Mn = np.array([[float(M[i,j]) for j in range(4)] for i in range(4)])
u, s, vh = np.linalg.svd(Mn)
print("singular values:", s)
ker = vh[-1]
print("numeric kernel vector:", ker, "ratios:", ker/ker[0])
# check candidate simple relations by hand: 2^a c0^b (1-c0)^c (1+c0)^d = +-1?
# relation lattice 1-dim; scale ker to integers
import flint
B = []
N = mp.mpf(10)**100
for i in range(4):
    r = [0]*8; r[i] = 1
    for j in range(4):
        r[4+j] = int(mp.nint(M[j,i]*N))
    B.append(r)
R = flint.fmpz_mat(B).lll().tolist()
for r in R[:3]:
    e = r[:4]
    res = max(abs(sum(e[i]*M[j,i] for i in range(4))) for j in range(4))
    print("LLL cand:", e, " resid:", mp.nstr(res, 3))
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=300)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
