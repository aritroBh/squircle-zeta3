t = r'''
import sympy as sp, mpmath as mp
import flint
c0 = sp.sin(sp.pi/5)
xx = sp.Symbol('xx')
mpoly = sp.minpoly(c0, xx)
roots = sp.real_roots(sp.Poly(mpoly, xx))
for prec in (60, 120, 240):
    mp.mp.dps = prec
    cvals = [mp.mpf(str(r.evalf(prec - 5))) for r in roots]
    M = mp.matrix([[mp.log(abs(v)) for v in (mp.mpf(2), cv, 1-cv, 1+cv)] for cv in cvals])
    print(f"prec {prec}: det = {mp.nstr(mp.det(M), 8)}")

# relation hunt at scale 10^250
mp.mp.dps = 280
cvals = [mp.mpf(str(r.evalf(270))) for r in roots]
M = mp.matrix([[mp.log(abs(v)) for v in (mp.mpf(2), cv, 1-cv, 1+cv)] for cv in cvals])
N = mp.mpf(10)**250
B = []
for i in range(4):
    r = [0]*8; r[i] = 1
    for j in range(4):
        r[4+j] = int(mp.nint(M[j,i]*N))
    B.append(r)
R = flint.fmpz_mat(B).lll().tolist()
for r in R[:4]:
    e = [int(x) for x in r[:4]]
    res = max(abs(mp.fsum(e[i]*M[j,i] for i in range(4))) for j in range(4))
    print("cand:", e, " resid:", mp.nstr(res, 3))
    if res < mp.mpf('1e-200'):
        g = [sp.Integer(2), c0, 1-c0, 1+c0]
        prod = sp.Integer(1)
        for gi, ei in zip(g, e):
            prod = sp.simplify(prod * gi**ei)
        print("   exact check prod =", prod, " => relation" if prod in (1,-1) else "   NOT exact")
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
