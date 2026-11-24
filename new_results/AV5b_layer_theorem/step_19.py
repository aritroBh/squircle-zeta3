
import mpmath as mp
mp.mp.dps = 80

P3_m3 = mp.mpf("-5.145455017135615837172")  # from step_18 (70 dps)
# recompute at higher dps to be safe
exec(open("step_18.py").read().replace("mp.mp.dps = 70","mp.mp.dps = 100").split("print(")[0])
P3_m3 = P_deriv(3, 3)
P2_m3 = P_deriv(3, 2)
P1_m3 = P_deriv(3, 1)
print("P3 m=3:", mp.nstr(P3_m3, 60))
print("P2 m=3:", mp.nstr(P2_m3, 60))

z3 = mp.zeta(3); pi2 = mp.pi**2; ln2 = mp.log(2)
basis = [z3, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
for scale, name in [(mp.sqrt(3), "sqrt3"), (mp.mpf(1), "1")]:
    x = P3_m3/scale
    # PSLQ
    M = [x] + basis
    rel = mp.pslq(M, tol=mp.mpf('1e-85'), maxcoeff=mp.mpf('1e15'), maxsteps=10000)
    print(f"scale {name}: PSLQ rel =", rel)
    if rel:
        resid = mp.fsum(r*v for r, v in zip(rel, M))
        print("   residual:", mp.nstr(resid, 6))
        print("   => P3 = %s*(%s)" % (name, str(rel)))
# also P2
x = P2_m3/mp.sqrt(3)
M = [x] + basis
rel = mp.pslq(M, tol=mp.mpf('1e-85'), maxcoeff=mp.mpf('1e15'))
print("P2/sqrt3 PSLQ:", rel)
if rel: print("   residual:", mp.nstr(mp.fsum(r*v for r,v in zip(rel,M)), 6))
