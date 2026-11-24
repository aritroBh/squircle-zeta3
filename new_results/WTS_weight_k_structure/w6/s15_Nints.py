"""w6 identification step 8: the 12 basic integrals
    N[(i,j),r] = int_b^a  dv/(v-r)  log^i t  log^j(1-t),   (i,j) in {(5,0),(4,1),(3,2)},  r in {0, a, -a, -b},
t=(v^2+2v-1)/(4v), a=1+sqrt2, b=sqrt2-1.  Same quadrature as s02_eval.py (tanh-sinh, split at v=1, integrating in the exact
distance to the endpoint).  Then W_6(k) = sum_{key,r} (rho0 + sqrt2 rho1)(k) * N  with the exact polynomials of s14.
usage: s15_Nints.py DPS ; writes Nints_DPS.json.   All printed numbers come from this script."""
import sys, os, json, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
DPS = int(sys.argv[1])
mp.mp.dps = DPS            # after imports
print("mp.mp.dps =", mp.mp.dps, flush=True)
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
POLE = {"0": mp.mpf(0), "a": aN, "-a": -aN, "-b": -bN}
def integ(i, j, r):
    def fl(s):
        v = bN + s; t = s * (v + aN) / (4 * v)
        return mp.log(t) ** i * mp.log1p(-t) ** j / (v - r)
    def fr(u):
        v = aN - u; w = u * (v + bN) / (4 * v)          # w = 1 - t
        return mp.log1p(-w) ** i * mp.log(w) ** j / ((-u) if r == aN else (v - r))
    v1, e1 = mp.quad(fl, [0, 1 - bN], error=True, maxdegree=12)
    v2, e2 = mp.quad(fr, [0, aN - 1], error=True, maxdegree=12)
    return v1 + v2, max(e1, e2)
out = {"dps": DPS, "N": {}}
t0 = time.time()
for (i, j) in ((5, 0), (4, 1), (3, 2)):
    for name, r in POLE.items():
        val, err = integ(i, j, r)
        out["N"][f"{i},{j}|{name}"] = mp.nstr(val, DPS)
        print(f"N[({i},{j}), {name}] = {mp.nstr(val, 30)}  err_est {mp.nstr(err, 3)} [{time.time()-t0:.0f}s]", flush=True)
json.dump(out, open(os.path.join(W, f"Nints_{DPS}.json"), "w"), indent=1)
