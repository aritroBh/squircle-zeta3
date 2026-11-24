import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J = json.load(open('out/d29/candidates_S7_B5.json'))
allc=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names)
for hmax,tag in ((2,'S7h2'),(4,'S7h4')):
    cands=[c for c in allc if c['h']<=hmax]
    K, comps, sym_dicts, rank = w5pipe.kernel_and_comps(cands, d, tag)
    json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank, hmax=hmax, saturated=False),
              open(f'out/d29/kernel5_{tag}.json','w'))
    print(f"saved kernel5_{tag}.json: N={len(cands)} kernel {len(K)} rank {rank}", flush=True)
    with open('progress.log','a') as f:
        f.write(f"step_15_d29{tag}_kernel: N={len(cands)} kernel {len(K)} rank {rank}\n")
