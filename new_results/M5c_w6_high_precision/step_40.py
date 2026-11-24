import datetime, pickle, math
work='<project>/research/08_theory/M5c_w6_high_precision'
def log(s):
    with open(work+'/progress.log','a') as f:
        f.write(f'{datetime.datetime.now().isoformat()} {s}\n')
log('leaves+atoms at 3500 digits done (99s); controls: Li6 dupl 4.4e-3611/2.8e-3611, Li32 vs M4c 9.5e-3542, N6 vs M5a <=4.6e-77')
# Two-precision leaf check: cached (ctx 12092 bits) vs fresh at 12700 bits
import sys
sys.dont_write_bytecode=True; sys.path.insert(0,work)
import mpmath as mp
from flint import ctx
from w6_values import leaf_value, leaf_recipes
recs=leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
if isinstance(saved,dict): saved=list(saved.values())
leaves={r: mp.mpf(v) for r,v in zip(recs,saved)}
ctx.prec=12700
mp.mp.dps=3700
maxdiff=mp.mpf(0); argmax=None
import time
t0=time.monotonic()
for rec in recs:
    v=mp.mpf(leaf_value(rec, 3620).mid().str(3680, radius=False))
    d=abs(v-leaves[rec])
    if d>maxdiff: maxdiff, argmax=d, rec
print('two-precision leaf check: recomputed at ctx=12700 bits vs cache(ctx=12092), %.1fs'%(time.monotonic()-t0))
print('max |diff| over %d leaves = %s at %s'%(len(recs), mp.nstr(maxdiff,6), argmax))
log(f'two-precision leaf check (12700 vs 12092 bits): max diff {mp.nstr(maxdiff,6)} at {argmax}')
