import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
from fractions import Fraction as Fr
mp.mp.dps = 1120
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2; pi4 = mp.pi**4
Li = lambda s, x: mp.polylog(s, x)
z3 = mp.zeta(3)
t0 = time.time()
Li4_args = {"Li4(-phi^3)": -phi**3, "Li4(-sqrt5 phi)": -r5*phi, "Li4(-2phi^2*phi)": -8-4*r5,
            "Li4(1/(2phi^2))": 1/(2*phi**2), "Li4(1/phi)": 1/phi, "Li4(1/phi^2)": 1/phi**2,
            "Li4(1/phi^3)": 1/phi**3, "Li4(1/2)": mp.mpf(1)/2, "Li4(2/phi^2)": 2/phi**2,
            "Li4(phi/2)": phi/2, "Li4(1/phi^6)": 1/phi**6, "Li4(-1/phi^3)": -1/phi**3}
li4 = {nm: Li(4, x) for nm, x in Li4_args.items()}
Li3_args = {"Li3(1/phi)": 1/phi, "Li3(1/phi^3)": 1/phi**3, "Li3(1/(2phi))": 1/(2*phi), "Li3(2/phi^2)": 2/phi**2}
li3 = {nm: Li(3, x) for nm, x in Li3_args.items()}
li2_13 = Li(2, 1/phi**3)
B4_NAMES = list(Li4_args) + ["L*zeta3", "lam*zeta3"]
for nm in Li3_args: B4_NAMES += ["L*"+nm, "lam*"+nm]
B4_NAMES += ["pi2*Li2(1/phi^3)", "Li2(1/phi^3)^2", "Li2(1/phi^3)*L^2", "Li2(1/phi^3)*L*lam", "Li2(1/phi^3)*lam^2",
             "pi4", "pi2*L^2", "pi2*L*lam", "pi2*lam^2", "L^4", "L^3*lam", "L^2*lam^2", "L*lam^3", "lam^4"]
B4 = list(li4.values()) + [L*z3, lam*z3]
for nm in Li3_args: B4 += [L*li3[nm], lam*li3[nm]]
B4 += [pi2*li2_13, li2_13**2, li2_13*L**2, li2_13*L*lam, li2_13*lam**2,
       pi4, pi2*L**2, pi2*L*lam, pi2*lam**2, L**4, L**3*lam, L**2*lam**2, L*lam**3, lam**4]
nb = len(B4)
print("W4 basis size:", nb, f"(eval {time.time()-t0:.0f}s)", "lattice dim", 2*nb+1)
import pickle
pickle.dump({"names": B4_NAMES, "vals": [mp.nstr(x, 1130) for x in B4]}, open(WD+"golden_B4.pkl", "wb"))
for tag in ("gol1", "gol3"):
    W4 = mp.mpf(open(WD+f"{tag}_W4_hpl_1150.txt").read())
    vals = [W4] + B4 + [r5*x for x in B4]
    cands = {}
    for dps in (700, 1000):
        t1=time.time()
        rows = L7.find_relation_rows(vals, dps)
        mp.mp.dps = 1120
        row = next(r for r in rows if r[0] != 0)
        cands[dps] = L7.normalize(row)
        resid = sum(c*v for c, v in zip(cands[dps], vals))
        print(f"{tag} dps {dps}: LLL {time.time()-t1:.0f}s, height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
    print(tag, "identical:", cands[700] == cands[1000])
    if cands[700] == cands[1000]:
        a = cands[1000]
        json.dump({"a0": a[0], "coeffs": a[1:], "basis": B4_NAMES}, open(WD+f"{tag}_W4_relation.json", "w"))
        print(f"  saved (a0 = {a[0]})")
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_43: golden W4 LLL attempt 2 (12 Li4 args, 73-dim).\n")
