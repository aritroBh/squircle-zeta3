"""Independent LLL identification of c_2, c_1, c_0 (python-flint LLL, gram='exact') against the 58-element Q-basis
{29 weight-4 atoms} x {1, sqrt2}, at several (data dps, scaling digits D) pairs.  Usage: nice -n 10 python id_final_lll.py <dps> <D> [targets...]
PRECISION: mp.mp.dps set after all imports and printed; a relation is accepted only if its residual (recomputed in mpmath at
mp.mp.dps) is below 10^-(D-10) and it is identical across the (dps, D) pairs (compared by id_compare.py)."""
import sys, json, time
from id_lib import *
dps = int(sys.argv[1]); D = int(sys.argv[2]); tg = sys.argv[3:] or ["c2", "c1", "c0", "delta"]
set_dps(dps); print("mp.mp.dps =", mp.mp.dps, " D =", D)
I = load_I(dps)
c = {f"c{j}": combo(I, f"c{j}") for j in range(3)}
c["delta"] = c["c0"] + 3*c["c1"] + 9*c["c2"] - combo(I, "k3")    # P(3) - W_4(3)
S = with_sqrt2(atoms()); labels = list(S.keys())
out = {}
for nm in tg:
    vals = [c[nm]] + list(S.values())
    t0 = time.time()
    res = find_relation(vals, D)
    co, s = res[0]
    co2, s2 = res[1]
    rel = {l: x for l, x in zip(["target"] + labels, co) if x}
    ok = abs(s) < mp.mpf(10)**(-(D - 10))
    print(f"{nm}: n={len(vals)} time {time.time()-t0:.1f}s  max|coef|={max(abs(x) for x in co)} (2nd row max|coef|={max(abs(x) for x in co2)})  residual={mp.nstr(s,4)}  accepted={ok}  mp.dps={mp.mp.dps}")
    print("   ", rel)
    out[nm] = {"relation": rel, "residual": mp.nstr(s, 5), "max_coef": max(abs(x) for x in co), "second_row_max_coef": max(abs(x) for x in co2), "accepted": ok}
fn = os.path.join(RES, f"id_lll_dps{dps}_D{D}.json")
prev = json.load(open(fn))["results"] if os.path.exists(fn) else {}      # merge, never clobber earlier targets
prev.update(out)
json.dump({"dps": dps, "D": D, "mp_dps": mp.mp.dps, "results": prev}, open(fn, "w"), indent=1)
