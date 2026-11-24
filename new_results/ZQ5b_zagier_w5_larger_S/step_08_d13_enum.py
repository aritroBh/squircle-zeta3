
import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
import os
os.makedirs('out/d13', exist_ok=True)
F = W.make_field(13, split=[3], inert=[2])
t0=time.time()
cands = W.enum_exceptional_fast(F, B=8)
data=[dict(x=[str(c['x'].a),str(c['x'].b)], vx=c['vx'], v1m=c['v1m'], h=c['h'], xf=float(c['x'].emb(40))) for c in cands]
json.dump(dict(D=13, B=8, split=[3], inert=[2], names=F['names'], disp=[str(W.fmt_qd(g)) for g in F['letters']],
               S="S={13,2,3} (13 ram, 2 inert, 3 split: (5+-sqrt13)/2)", cands=data), open('out/d13/candidates.json','w'), indent=1)
print(f"D=13: saved {len(data)} exceptional S-units (B=8) ({time.time()-t0:.0f}s)")
import collections
print("h dist:", dict(sorted(collections.Counter(c['h'] for c in cands).items())))
