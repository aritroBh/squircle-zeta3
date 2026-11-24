import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
import importlib, cor9_rank; importlib.reload(cor9_rank)
from cor9_rank import true_rank
survivors = pickle.load(open(WD+"census_filter_survivors.pkl","rb"))
t0 = time.time()
rank2 = []
fails = []
ntot = 0
for D in sorted(survivors):
    for a, b, e in survivors[D]:
        ntot += 1
        r, M, pids, rel = true_rank(D, a, b, e)
        if r is None:
            fails.append((D, a, b, e))
        elif r <= 2:
            rank2.append((D, a, b, e, r))
    print(f"D={D} done ({time.time()-t0:.0f}s), rank2: {len(rank2)}, fails: {len(fails)}", flush=True)
print("\nRANK-2 (complete exact census, D<=200, heights<=40):", len(rank2))
import mpmath as mp
mp.mp.dps = 50
for D, a, b, e, r in rank2:
    c1 = (mp.mpf(a) + mp.mpf(b)*mp.sqrt(D))/e
    print(f"  Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}s{D})/{e} = {mp.nstr(c1,20)}  alpha = {mp.nstr(2*mp.acos(c1),20)}  rank {r}")
print("unit-check failures:", fails)
json.dump(rank2, open(WD+"census_rank2_complete.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_62: complete exact census via true_rank; {len(rank2)} rank-2 points.\n")
