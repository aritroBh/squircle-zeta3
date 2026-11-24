# Identification step 4: c_3 (leading k-coefficient of W_5(k)) in the Q-span of weight-5 alternating-MZV constants.
# Values: results/cm_D1500_m10.json (1500 digits stored, 1350 verified by m=9 vs m=10 step-halving).
# PSLQ at two (several) precisions; relation accepted only if identical at all of them.
import sys, os, json, time
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 1450                      # AFTER imports; reset per run below
R = json.load(open("results/cm_D1500_m10.json"))
mp.mp.dps = 1450
c3 = mp.mpf(R["c"][3])
lam = mp.log(2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5)
Li4h = mp.polylog(4, mp.mpf(1)/2); Li5h = mp.polylog(5, mp.mpf(1)/2)
names8 = ["zeta5", "pi2 zeta3", "zeta3 lam2", "pi4 lam", "pi2 lam3", "lam5", "Li4(1/2) lam", "Li5(1/2)"]
vals8 = [z5, pi**2*z3, z3*lam**2, pi**4*lam, pi**2*lam**3, lam**5, Li4h*lam, Li5h]
log = {}
for D in (300, 285, 600, 585, 1300, 1285):
    mp.mp.dps = D + 20
    # recompute the basis at this working precision from scratch (constants are cheap)
    lam = mp.log(2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5)
    Li4h = mp.polylog(4, mp.mpf(1)/2); Li5h = mp.polylog(5, mp.mpf(1)/2)
    vals8 = [z5, pi**2*z3, z3*lam**2, pi**4*lam, pi**2*lam**3, lam**5, Li4h*lam, Li5h]
    x = mp.mpf(R["c"][3])
    t0 = time.time()
    rel = mp.pslq([x] + vals8, tol=mp.mpf(10)**(-D), maxcoeff=10**12, maxsteps=10**6)
    res = sum(c*v for c, v in zip(rel, [x] + vals8)) if rel else None
    print(f"D={D} mp.mp.dps={mp.mp.dps}: rel={rel}  residual={mp.nstr(res, 5) if res is not None else None}  ({time.time()-t0:.1f}s)", flush=True)
    log[D] = rel
json.dump({str(k): v for k, v in log.items()}, open("results/id4_c3_pslq.json", "w"))
