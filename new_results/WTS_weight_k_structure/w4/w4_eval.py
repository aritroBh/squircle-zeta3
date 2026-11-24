"""Usage: nice -n 10 python w4_eval.py <weight> <k> <dps>
Evaluates W_weight(k) = sum_{i+j=weight-1} int_b^a K_ij log^i t log^j(1-t) dv by tanh-sinh (split at v=1).
Saves results/w<weight>_k<k>_dps<dps>.json (resumable: skips if present)."""
import sys, os, json, time
from w4_common import *
w, k, dps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
fn = os.path.join(RES, f"w{w}_k{k}_dps{dps}.json")
if os.path.exists(fn):
    print("exists", fn); sys.exit(0)
set_dps(dps)                                   # after all imports
print("mp.mp.dps =", mp.mp.dps, flush=True)
d = w - 1
pcs = pieces_mp(k, d)
t0 = time.time()
total = mp.mpf(0); errtot = mp.mpf(0); per = {}
for (i, j, poles) in pcs:
    val, err = integrate([(i, j, poles)])
    per[f"{i},{j}"] = {"value": mp.nstr(val, dps + 5), "quad_err": mp.nstr(err, 5)}
    total += val; errtot += err
    print(f"  w{w} k={k} key({i},{j}) = {mp.nstr(val, 30)}  quad_err {mp.nstr(err, 5)}  dps {mp.mp.dps}", flush=True)
json.dump({"weight": w, "k": k, "dps": dps, "mp_dps_at_end": mp.mp.dps, "value": mp.nstr(total, dps + 5),
           "quad_err_sum": mp.nstr(errtot, 5), "per_key": per, "seconds": time.time() - t0},
          open(fn, "w"), indent=1)
print(f"W_{w}({k}) = {mp.nstr(total, 40)}   [{time.time()-t0:.1f}s, dps {mp.mp.dps}]")
