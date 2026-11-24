import math, pickle, time
import sys
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3700
from flint import arb, acb, ctx
from w6_values import numeric, leaf_recipes
from w4_symbol import Q, R, A, B
from fractions import Fraction as F
ONE=Q(1)

def mpl_general2(s,t,x,y,digits):
    """Li_{s,t}(x,y) with N from the TRUE decay rate rx*max(1,ry) and rigorous tail."""
    xx, yy = numeric(x), numeric(y)
    rx, ry = abs(xx), abs(yy)
    assert 0 < rx < 1 and rx*ry < 1
    rate = rx*(ry if ry>1 else 1)
    N = int((digits+40)*math.log(10)/-math.log(float(rate))) + 60
    total, H, xp, yp = arb(0), arb(0), arb(1), arb(1)
    for n in range(1, N+1):
        xp *= xx
        total += xp*H/(n**s)
        yp *= yy
        H += yp/(n**t)
    if ry > 1:
        tail = (rx*ry)**(N+1)/((ry-1)*(N+1)**s*(1-rx*ry))
    else:
        Hb = arb(2).zeta() if t==2 else arb(3).zeta()  # |H(n)| <= zeta(t)
        tail = Hb*rx**(N+1)/((N+1)**s*(1-rx))
    return total + arb(0, tail.abs_upper()), N

ctx.prec = 12700
mpls = [r for r in leaf_recipes() if r[0]=='mpl']
saved = pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
recs = leaf_recipes()
old = {r: mp.mpf(v) for r,v in zip(recs, saved)}
t0=time.monotonic()
for rec in mpls:
    _,s,t,x,y = rec
    v, N = mpl_general2(s,t,x,y,3560)
    nv = mp.mpf(v.mid().str(3640, radius=False))
    rad = mp.mpf(v.rad().str(6))
    print(rec, 'N=%d'%N, 'rad=%s'%v.rad().str(6), 'old-new=%s'%mp.nstr(abs(nv-old[rec]),5), '(%.1fs)'%(time.monotonic()-t0))
