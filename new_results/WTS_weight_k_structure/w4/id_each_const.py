"""Optional: LLL identification of each of the 8 fixed constants J_r, K_r against the 58-dim basis (dps 600, D 585)."""
import json, time
from id_lib import *
dps, D = 600, 585
set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
I = load_I(dps); S = with_sqrt2(atoms()); labels = list(S.keys())
out = {}
for n in sorted(I):
    t0 = time.time()
    res = find_relation([I[n]] + list(S.values()), D)
    co, s = res[0]; co2, _ = res[1]
    ok = abs(s) < mp.mpf(10)**(-(D - 10))
    print(f"{n}: max|coef|={max(abs(x) for x in co)}  2nd row {max(abs(x) for x in co2)}  residual {mp.nstr(s,4)}  target coef {co[0]}  accepted={ok}  nterms={sum(1 for x in co if x)}  [{time.time()-t0:.0f}s] dps {mp.mp.dps}", flush=True)
    out[n] = {"relation": {l: x for l, x in zip(["target"] + labels, co) if x}, "max_coef": max(abs(x) for x in co), "second_row": max(abs(x) for x in co2), "residual": mp.nstr(s, 5), "accepted": ok}
json.dump(out, open(os.path.join(RES, "id_each_const_lll.json"), "w"), indent=1)
