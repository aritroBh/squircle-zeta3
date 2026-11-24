import sys, time, pickle, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import trunc_engine as te
from fastred import QSI, ZERO
def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def build_tower_min(k, D, rho_k, nmin):
    tower = {key: list(vals) for key, vals in rho_k.items()}
    for n in range(k-1, nmin-1, -1):
        dmax_n = min(D, 2*n-4)
        for key in te.keys_of_degree(n):
            i, j = key
            def U(kk, dp, _t=tower):
                v = _t.get(kk)
                if v is None: return ZERO
                return v[dp] if dp < len(v) else ZERO
            tower[key] = [te.Lambda_ij(i, j, delta, U) for delta in range(dmax_n+1)]
    return tower
t0=time.time()
S2 = te.stirling_table(32)
TOWERS = {}
DEGMAX = 18
for k in range(7, 31):
    rho_k = te.rho_level(k, 10, S2)
    tw = build_tower_min(k, 10, rho_k, 7)
    for n in range(7, DEGMAX+1):
        if n > k: continue
        for key in te.keys_of_degree(n):
            TOWERS.setdefault(k, {})[key] = {delta: q4(tw[key][delta]) for delta in range(0, min(10, 2*n-4)+1)}
pickle.dump(TOWERS, open(work+'/data_towers_D10_K30_deg18.pkl','wb'))
# count cells
totalcells = sum(1 for k in TOWERS for key in TOWERS[k] for delta in TOWERS[k][key])
n9 = sum(1 for k in TOWERS for key in TOWERS[k] if 9 in TOWERS[k][key])
n10 = sum(1 for k in TOWERS for key in TOWERS[k] if 10 in TOWERS[k][key])
print(f"extended towers k=7..30 deg 7..18 D=10 ({round(time.time()-t0,1)}s); depth9 {n9}, depth10 {n10}, total {totalcells} cells")
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} stepB3: extended towers k=7..30 deg 7..18; {totalcells} cells\n")
