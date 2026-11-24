"""w6 identification step 9: reconstruct c_m (m=0..4) = sum rho^{(m)}_{key,r,comp} * N_{key,r} * sqrt2^comp from the exact rho (s14) and the
integrals N (s15), and compare with the fitted c_m (s08).  Controls the whole new pipeline.  Printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from fractions import Fraction as Fr
rho = json.load(open(os.path.join(W, "rho_exact.json")))
for dps, run in ((250, 250), (400, 400)):
    mp.mp.dps = dps
    r2 = mp.sqrt(2)
    N = {k: mp.mpf(v) for k, v in json.load(open(os.path.join(W, f"Nints_{dps}.json")))["N"].items()}
    for m in range(5):
        tot = mp.mpf(0)
        for name, co in rho.items():
            key, let, comp = name.split("|")
            f = Fr(co[m])
            if f == 0: continue
            tot += mp.mpf(f.numerator) / mp.mpf(f.denominator) * N[f"{key}|{let}"] * (r2 if comp == "1" else 1)
        print(f"dps {mp.mp.dps}: c_{m}(from N) - c_{m}(fit, run {run}) = {mp.nstr(tot - mp.mpf(coef(run, m)), 5)}")
