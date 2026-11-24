import time, json, pickle
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
from fractions import Fraction as Fr
mp.mp.dps = 1120
r5 = mp.sqrt(5)
d = pickle.load(open(WD+"golden_B4.pkl", "rb"))
B4_NAMES = d["names"]; B4 = [mp.mpf(v) for v in d["vals"]]
keep = [i for i, nm in enumerate(B4_NAMES) if nm not in
        ('Li4(-phi^3)', 'Li4(-2phi^2*phi)', 'Li4(1/phi)', 'Li4(1/phi^3)')]
B4r_NAMES = [B4_NAMES[i] for i in keep]; B4r = [B4[i] for i in keep]
nb = len(B4r)
vals = B4r + [r5*x for x in B4r]
rows = L7.find_relation_rows(vals, 900)
mp.mp.dps = 1120
pure = sum(1 for r in rows if abs(sum(c*v for c, v in zip(L7.normalize(r), vals))) < mp.mpf('1e-1000'))
print("pure-basis relations in reduced 32-el basis:", pure)
for tag in ("gol1", "gol3"):
    W4 = mp.mpf(open(WD+f"{tag}_W4_hpl_1150.txt").read())
    vals2 = [W4] + B4r + [r5*x for x in B4r]
    cands = {}
    for dps in (700, 1000):
        t1=time.time()
        rows = L7.find_relation_rows(vals2, dps)
        mp.mp.dps = 1120
        row = next(r for r in rows if r[0] != 0)
        cands[dps] = L7.normalize(row)
        resid = sum(c*v for c, v in zip(cands[dps], vals2))
        print(f"{tag} dps {dps}: LLL {time.time()-t1:.0f}s, height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
    print(tag, "identical:", cands[700] == cands[1000])
    if cands[700] == cands[1000]:
        a = cands[1000]
        print(f"W4({tag}), a0 = {a[0]}:")
        for i in range(1, nb+1):
            p, q = -a[i], -a[i+nb]
            if p or q:
                print(f"  ({Fr(p, a[0])} {'+' if q>=0 else ''}{Fr(q, a[0])}*sqrt5) * {B4r_NAMES[i-1]}")
        json.dump({"a0": a[0], "coeffs": a[1:], "basis": B4r_NAMES}, open(WD+f"{tag}_W4_relation.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_45: golden W4 LLL with independent 32-el basis.\n")
