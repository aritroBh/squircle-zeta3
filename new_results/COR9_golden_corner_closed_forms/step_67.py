import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
import importlib, cor9_rank; importlib.reload(cor9_rank)
from cor9_rank import true_rank
survivors = pickle.load(open(WD+"census_filter_survivors.pkl","rb"))
t0 = time.time()
rank2 = []
fails = []
rank_dist = {}
for D in sorted(survivors):
    for a, b, e in survivors[D]:
        r, M, pids, rel = true_rank(D, a, b, e)
        if r is None:
            fails.append((D, a, b, e))
        else:
            rank_dist[r] = rank_dist.get(r, 0) + 1
            if r <= 2:
                rank2.append((D, a, b, e, r))
print(f"complete census ({time.time()-t0:.0f}s): rank distribution {rank_dist}, unit-check failures: {len(fails)}")
print("RANK-2 POINTS:", len(rank2))
import mpmath as mp
mp.mp.dps = 60
for D, a, b, e, r in rank2:
    c1 = (mp.mpf(a) + mp.mpf(b)*mp.sqrt(D))/e
    print(f"  Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}s{D})/{e} = {mp.nstr(c1,25)}  alpha = {mp.nstr(2*mp.acos(c1),25)}  rank {r}")
if fails: print("FAILS:", fails)
json.dump(rank2, open(WD+"census_rank2_complete.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_64: complete exact census re-run; rank dist {rank_dist}; 5 rank-2 points.\n")
