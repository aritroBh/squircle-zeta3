"""verify step: own W_{d+1}(k) values (own Gauss-Legendre quadrature of the exact finals_k.json).
usage: verify_vals.py DPS OUTDIR D K [K ...]      (D = d, weight = d+1)
Writes OUTDIR/own_k{K}_d{D}_{DPS}.json (resumable: skips existing).  mp.mp.dps = DPS + 10 is set after imports and printed.
"""
import sys, os, json, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp

DPS = int(sys.argv[1]); OUT = sys.argv[2]; D = int(sys.argv[3]); KS = [int(x) for x in sys.argv[4:]]
os.makedirs(OUT, exist_ok=True)
WD = DPS + 10
VL.set_dps(WD)
print(f"mp.mp.dps = {mp.mp.dps} (target {DPS} digits)", flush=True)
for k in KS:
    fn = os.path.join(OUT, f"own_k{k}_d{D}_{DPS}.json")
    if os.path.exists(fn):
        print(f"k={k} d={D}: exists, skipping", flush=True)
        continue
    t0 = time.time()
    fin = VL.load_finals(os.path.join(W, f"finals_{k}.json"))
    val, pcs = VL.weight_part(fin, D, WD)
    rec = {"k": k, "d": D, "dps_target": DPS, "mp_mp_dps": mp.mp.dps, "value": mp.nstr(val, DPS + 5),
           "pieces": {f"{i},{j}": mp.nstr(v, DPS + 5) for (i, j), v in pcs.items()}, "seconds": time.time() - t0}
    tmp = fn + ".tmp"
    json.dump(rec, open(tmp, "w"), indent=1)
    os.replace(tmp, fn)
    print(f"k={k} d={D}: W_{D+1} = {mp.nstr(val, 50)}   [{time.time()-t0:.0f}s] mp.mp.dps={mp.mp.dps}", flush=True)
