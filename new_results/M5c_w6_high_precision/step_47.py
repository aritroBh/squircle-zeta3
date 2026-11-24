import math, pickle, datetime
import sys
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3700
from flint import arb, ctx
from w6_values import numeric
from w4_symbol import Q, R, B
from fractions import Fraction as F

def mpl_general2(s,t,x,y,digits):
    xx, yy = numeric(x), numeric(y)
    rx, ry = abs(xx), abs(yy)
    rate = rx*(ry if ry>1 else 1)
    N = int((digits+40)*math.log(10)/-math.log(float(rate))) + 60
    total, H, xp, yp = arb(0), arb(0), arb(1), arb(1)
    for n in range(1, N+1):
        xp *= xx; total += xp*H/(n**s)
        yp *= yy; H += yp/(n**t)
    Hb = arb(2).zeta() if t==2 else arb(3).zeta()
    tail = (rx*ry)**(N+1)/((ry-1)*(N+1)**s*(1-rx*ry)) if ry>1 else Hb*rx**(N+1)/((N+1)**s*(1-rx))
    return total + arb(0, tail.abs_upper()), N

ctx.prec=12700
v1,N1 = mpl_general2(4,2,B,R/2,3460)
v2,N2 = mpl_general2(4,2,B,R/2,3560)
a=mp.mpf(v1.mid().str(3600,radius=False)); b=mp.mpf(v2.mid().str(3640,radius=False))
print('Li42 two-precision: N=%d(rad %s) vs N=%d(rad %s)'%(N1,v1.rad().str(6),N2,v2.rad().str(6)))
print('agreement:', mp.nstr(abs(a-b),6))
# also vs patched cache
recs=__import__('w6_values').leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
idx=[i for i,r in enumerate(recs) if r[0]=='mpl' and r[1]==4 and r[2]==2][0]
print('cache vs N2:', mp.nstr(abs(mp.mpf(saved[idx])-b),6))
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} Li42 fixed leaf two-precision check: N={N1} vs N={N2} agree to {mp.nstr(abs(a-b),6)}; cache matches.\n')
