import sys, pickle, time
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3700
from flint import ctx
from w6_values import leaf_value, leaf_recipes
recs=leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
leaves={r: mp.mpf(v) for r,v in zip(recs,saved)}
ctx.prec=12700
t0=time.monotonic()
worst=[]
for rec in recs:
    v=mp.mpf(leaf_value(rec, 3580).mid().str(3640, radius=False))
    d=abs(v-leaves[rec])
    worst.append((d,rec))
worst.sort(reverse=True)
for d,rec in worst[:12]:
    print(mp.nstr(d,4), rec)
print('(%.0fs, %d leaves)'%(time.monotonic()-t0, len(recs)))
