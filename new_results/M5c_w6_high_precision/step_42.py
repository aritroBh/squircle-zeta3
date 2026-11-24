import datetime, pickle, math, time
work='<project>/research/08_theory/M5c_w6_high_precision'
import sys
sys.dont_write_bytecode=True; sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3700   # BEFORE loading pickle!
from flint import ctx
from w6_values import leaf_value, leaf_recipes
recs=leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
if isinstance(saved,dict): saved=list(saved.values())
leaves={r: mp.mpf(v) for r,v in zip(recs,saved)}
ctx.prec=12700
maxdiff=mp.mpf(0); argmax=None
t0=time.monotonic()
for rec in recs:
    v=mp.mpf(leaf_value(rec, 3580).mid().str(3640, radius=False))
    d=abs(v-leaves[rec])
    if d>maxdiff: maxdiff, argmax=d, rec
print('two-precision leaf check: recomputed at ctx=12700 bits vs cache(ctx=12092), %.1fs'%(time.monotonic()-t0))
print('max |diff| over %d leaves = %s at %s'%(len(recs), mp.nstr(maxdiff,6), argmax))
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} two-precision leaf check FIXED (12700 vs 12092 bits): max diff {mp.nstr(maxdiff,6)} at {argmax}\n')
