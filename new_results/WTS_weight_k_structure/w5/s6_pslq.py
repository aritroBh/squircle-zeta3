# Step 6: closed form of the leading coefficient c_3 (and exploratory PSLQ for the others).
import sys, os, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 140                     # after imports
print("mp.mp.dps =", mp.mp.dps)
R = json.load(open("s5_results.json"))
cm = [mp.mpf(x) for x in R["cm_direct"]]          # 120 stored digits
delta = mp.mpf(R["defect_k4"])
sq2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+sq2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5)
# closed form of c_3 from PSLQ relation [-3456, 270, 30, 9, 180, 20, 12]
c3 = (270*z5 + 30*pi**2*z3 + 9*pi**4*lam + 180*z3*lam**2 + 20*pi**2*lam**3 + 12*lam**5) / 3456
print("c_3 closed form - c_3 (direct kernel integral) =", mp.nstr(c3 - cm[3], 5), " [stored digits: 135; mp.mp.dps =", mp.mp.dps, "]")
# exploratory
B1n = ["z5", "pi2 z3", "pi4 lam", "z3 lam2", "pi2 lam3", "lam5"]
B1 = [z5, pi**2*z3, pi**4*lam, z3*lam**2, pi**2*lam**3, lam**5]
B2n = B1n + ["pi4 L", "z3 lam L", "z3 L2", "pi2 lam2 L", "pi2 lam L2", "pi2 L3", "lam4 L", "lam3 L2", "lam2 L3", "lam L4", "L5"]
B2 = B1 + [pi**4*L, z3*lam*L, z3*L**2, pi**2*lam**2*L, pi**2*lam*L**2, pi**2*L**3, lam**4*L, lam**3*L**2, lam**2*L**3, lam*L**4, L**5]
def tryit(label, x, names, basis, tol=mp.mpf(10)**-105, maxc=10**6):
    r = mp.pslq([x] + basis, tol=tol, maxcoeff=maxc, maxsteps=10**6)
    print(f"  {label}: relation = {r}")
    return r
out = {}
P = lambda s: cm[0] + cm[1]*s + cm[2]*s**2 + cm[3]*s**3
tests = {"c2": cm[2], "c1": cm[1], "c0": cm[0], "delta_k4_defect": delta}
for s in (3, 4, 5, 6): tests[f"c1+2*{s}*c2 (slope at k={s})"] = cm[1] + 2*s*cm[2] + 3*s*s*cm[3]
for s in (3, 4, 5): tests[f"curv at k={s}: c2+3*{s}*c3"] = cm[2] + 3*s*cm[3]
for name, x in tests.items():
    print(name, mp.nstr(x, 30))
    for bn, bs, bb in (("B1(Q-span of z5,pi2z3,pi4lam,z3lam2,pi2lam3,lam5)", B1n, B1), ("B1+sqrt2*B1", B1n, B1 + [sq2*t for t in B1]), ("B2(17, with L)", B2n, B2)):
        out[f"{name}|{bn}"] = tryit(bn, x, bs, bb)
json.dump(dict(c3_closed_form="(270 zeta(5) + 30 pi^2 zeta(3) + 9 pi^4 log2 + 180 zeta(3) log^2 2 + 20 pi^2 log^3 2 + 12 log^5 2)/3456",
               c3_resid=mp.nstr(c3 - cm[3], 5), exploratory={k: v for k, v in out.items()}), open("s6_results.json", "w"), indent=1, default=str)
