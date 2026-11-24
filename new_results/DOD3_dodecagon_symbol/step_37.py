import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 200
import sympy as sp
from dod3_R import *
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
g = pickle.load(open(os.path.join(d2,'gdswaps620.pkl'),'rb'))
SwG, SwD = g['SwG'], g['SwD']

Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def CreN(al,be,u):
    pts = [mp.mpf(0)]
    for kap in (al,be):
        if kap > 0:
            s0 = 1/kap
            if 0 < s0 < u: pts.append(s0)
    pts.append(u); pts = sorted(set(pts))
    tot = mp.mpf(0)
    for a_,b_ in zip(pts[:-1], pts[1:]):
        tot += mp.quad(lambda s: mp.log(abs(1-al*s))*mp.log(abs(1-be*s))/s, [a_,b_])
    return tot
def SwremN(ka,kb):
    return -mp.log(abs(1-ka))*(Li2R(kb)-Li2R(kb/ka)) + Li3R(kb)

for (i,j,k) in [(6,1,5),(1,2,5),(1,2,6),(2,5,6)]:
    xi = Xs[i]; a = pdiv(Xs[j],xi); b = pdiv(Xs[k],xi)
    ka = pdiv(a, psub(a,ONE)); kb = pdiv(b, psub(b,ONE))
    lo = psub(ONE,xi); da, db = pmul(ka,lo), pmul(kb,lo)
    kan,kbn,dan,dbn = pfloat(ka),pfloat(kb),pfloat(da),pfloat(db)
    b1 = (CreN(kan,kbn,mp.mpf(1)) - SwG[(i,j,k)]) - SwremN(kan,kbn)
    b2 = (CreN(dan,dbn,mp.mpf(1)) - SwD[(i,j,k)]) - SwremN(dan,dbn)
    print((i,j,k), " bracketG =", mp.nstr(b1,4), " bracketD =", mp.nstr(b2,4), " diff =", mp.nstr(b1-b2,4))
    print("     ka,kb =",mp.nstr(kan,6),mp.nstr(kbn,6)," da,db =",mp.nstr(dan,6),mp.nstr(dbn,6))
