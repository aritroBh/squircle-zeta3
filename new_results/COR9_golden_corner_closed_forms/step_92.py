import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
import importlib, cor9_rank; importlib.reload(cor9_rank)
from cor9_rank import true_rank, rank_int
survivors = pickle.load(open(WD+"census_filter_survivors_v2.pkl","rb"))
t0 = time.time()
rank2 = []
fails = []
rank_dist = {}
decomp = {}
for D in sorted(survivors):
    for a, b, e in survivors[D]:
        r, M, pids, rel = true_rank(D, a, b, e)
        if r is None:
            fails.append((D, a, b, e)); continue
        rank_dist[r] = rank_dist.get(r, 0) + 1
        rv = rank_int(M)
        oddcols = [j for j, pid in enumerate(pids) if pid[0] != 2]
        rvodd = rank_int([[row[j] for j in oddcols] for row in M])
        key = (r, rv, rvodd, r - rv)
        decomp[key] = decomp.get(key, 0) + 1
        if r <= 2:
            rank2.append((D, a, b, e, r))
print(f"corrected census ({time.time()-t0:.0f}s): rank distribution {rank_dist}, unit-check failures: {len(fails)}")
print("RANK-2 POINTS:", len(rank2))
import mpmath as mp
mp.mp.dps = 60
for D, a, b, e, r in rank2:
    c1 = (mp.mpf(a) + mp.mpf(b)*mp.sqrt(D))/e
    print(f"  Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}s{D})/{e} = {mp.nstr(c1,25)}  alpha = {mp.nstr(2*mp.acos(c1),25)}  rank {r}")
print("(total rank, rank V, rank V_odd, delta) -> count:")
for k in sorted(decomp): print("  ", k, "->", decomp[k])
if fails: print("FAILS:", fails)
json.dump(rank2, open(WD+"census_rank2_complete.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_90: corrected census on 862 survivors; rank dist {rank_dist}; 5 rank-2 points.\n")
