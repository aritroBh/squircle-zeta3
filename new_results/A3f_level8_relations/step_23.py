import sys, time
sys.dont_write_bytecode=True
sys.path.insert(0,'Q')
from flint import ctx, acb, arb
import check_conversion as CC
ctx.prec = 1100
ctx.cap = 2048
CC.N = 2048
CC.roots = [acb(0, arb.pi()*j/4).exp() for j in range(8)]
CC.roots[0], CC.roots[2], CC.roots[4], CC.roots[6] = acb(1), acb(0,1), acb(-1), acb(0,-1)
from check_conversion import h
import mpmath as mp
mp.mp.dps=335
t0=time.monotonic()
v=h((7,7))
t1=time.monotonic()
v2=h((7,0))
t2=time.monotonic()
print('one word time', t1-t0, 'second', t2-t1)
print('radius:', mp.nstr(mp.mpf(v.real.rad().str(30,radius=False,more=True)),8))
# check known relation residual directly
import pickle
rows = pickle.load(open('audit_weight2.pickle','rb'))['rows']
from check_conversion import evaluate, approx
worst=mp.mpf(0)
for row in rows:
    s = abs(approx(evaluate(dict(row))))
    worst=max(worst,s)
print('known rows worst residual (direct evaluate):', mp.nstr(worst,8), 'total time', time.monotonic()-t0)
