
import mpmath as mp
mp.mp.dps = 100
exec(open("step_18.py").read().replace("mp.mp.dps = 70","mp.mp.dps = 100").split("print(")[0])
P3_m3 = P_deriv(3, 3); P2_m3 = P_deriv(3, 2)
z3 = mp.zeta(3); pi2 = mp.pi**2; ln2 = mp.log(2)
basis = [z3, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
names = ["z3","pi2ln2","ln2^3","pi2","ln2^2","ln2","1"]
for scale, sname in [(mp.sqrt(3), "sqrt3"), (mp.mpf(1), "1"), (mp.mpf(3), "3")]:
    for val, vname in [(P3_m3,"P3"), (P2_m3,"P2")]:
        M = [val/scale] + basis
        rel = mp.pslq(M, tol=mp.mpf('1e-90'), maxcoeff=10**15, maxsteps=10000)
        if rel:
            resid = mp.fsum(mp.mpf(r)*v for r, v in zip(rel, M))
            print(f"{vname}/{sname}: rel={rel} resid={mp.nstr(resid,4)}")
        else:
            print(f"{vname}/{sname}: no relation")
