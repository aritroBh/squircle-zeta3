
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr, pslq
mp.dps = 110
exec(open('av2core.py').read())
z2, z3 = zeta(2), zeta(3)

print("=== B2(p) at p = 1/(1-tau), tau grid ===")
taus = [-0.1, -0.05, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.05, 0.1]
vals = []
for tau in taus:
    p = 1/(1-tau)
    v, _ = B2_series(p)
    vals.append((tau, v))
    print("tau=%6.3f p=%8.5f B2=%s" % (tau, p, nstr(v, 40)))
# fit B2 = b0 + d1 t + d2 t^2 + d3 t^3 + d4 t^4
A = []; b = []
for (tau, v) in vals:
    A.append([mpf(1), tau, tau**2, tau**3, tau**4])
    b.append(v)
Am = matrix(A); bm = matrix(b)
sol = lu_solve(Am.T*Am, Am.T*bm)
res = Am*sol-bm
print("fit residual:", nstr(max(fabs(x) for x in res),5))
b0, d1, d2, d3, d4 = sol
for nm, s in zip(['b0','d1','d2','d3','d4'], sol):
    print(nm, "=", nstr(s, 40))
print("b0 vs 2 zeta(3)/3:", nstr(b0-2*z3/3, 10))
def hunt(name, val, basis, tol=mpf('1e-40'), maxcoeff=10**8):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if rel is None or rel[0] == 0:
        print("%s: (none)" % name); return
    resid = fabs(sum(mpf(rel[i])*([val]+basis)[i] for i in range(len(rel))))
    scale = max(fabs(mpf(rel[i])*([val]+basis)[i]) for i in range(len(rel)))
    tag = "OK" if resid/scale < mpf('1e-35') else "WEAK"
    print("%s [%s]: %s  resid=%s" % (name, tag, rel, nstr(resid/scale, 5)))
hunt("d1", d1, [mpf(1), z2, z3, z2*z3, log(2), z3*log(2)])
hunt("d2", d2, [mpf(1), z2, z3, z2*z3, log(2), z3*log(2)])
hunt("d1b", d1, [mpf(1), z2, z3, zeta(4), z2**2, log(2), log(3)])
hunt("d2b", d2, [mpf(1), z2, z3, zeta(4), z2**2, log(2), log(3)])
