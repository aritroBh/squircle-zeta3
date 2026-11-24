"""The 8 fixed weight-4 constants I_{(i,j),r} = int_b^a dv/(v-r) log^i t log^j(1-t), (i,j) in {(3,0),(2,1)}, r in {-a,0,-b,a}.
Usage: nice -n 10 python id_const.py <dps>.  Saves results/id_const_dps<dps>.json (resumable).
PRECISION: mp.mp.dps set after all imports and printed."""
import sys, json, time
from w4_common import *
dps = int(sys.argv[1])
fn = os.path.join(RES, f"id_const_dps{dps}.json")
if os.path.exists(fn): print("exists", fn); sys.exit(0)
set_dps(dps); print("mp.mp.dps =", mp.mp.dps, flush=True)
a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
poles = {"-a": -a, "0": mp.mpf(0), "-b": -b, "a": a}
out = {}; t0 = time.time()
for (i, j) in ((2, 1), (3, 0)):
    for name, r in poles.items():
        f = make_integrand([(i, j, [(r, mp.mpf(1))])])
        val, err = mp.quad(f, [b, 1, a], error=True)
        out[f"{i}{j}|{name}"] = mp.nstr(val, dps + 5)
        print(f"  I[{i},{j}; r={name}] = {mp.nstr(val, 25)}  quad_err {mp.nstr(err, 3)}  [{time.time()-t0:.0f}s] dps {mp.mp.dps}", flush=True)
json.dump({"dps": dps, "mp_dps_at_end": mp.mp.dps, "I": out}, open(fn, "w"), indent=1)
