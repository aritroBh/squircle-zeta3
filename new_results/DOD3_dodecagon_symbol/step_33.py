import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 200
import sympy as sp
from dod3_R import *

Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def SwremN(ka,kb):
    return -mp.log(abs(1-ka))*(Li2R(kb)-Li2R(kb/ka)) + Li3R(kb)

vals200 = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals200)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,300))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,300)))))
    return mp.mpf(sp.N(e, 200))

for (i,j,k) in [(6,1,5),(1,2,5),(1,2,6),(2,5,6)]:
    xi = Xs[i]; a = pdiv(Xs[j],xi); b = pdiv(Xs[k],xi)
    ka = pdiv(a, psub(a,ONE)); kb = pdiv(b, psub(b,ONE))
    lo = psub(ONE,xi)
    da, db = pmul(ka,lo), pmul(kb,lo)
    s1s, s2s = Swrem(ka,kb), Swrem(da,db)
    s1n = SwremN(pfloat(ka),pfloat(kb)); s2n = SwremN(pfloat(da),pfloat(db))
    print((i,j,k))
    print("   Swrem(ka,kb): sym %s num %s"%(mp.nstr(ev(s1s),14), mp.nstr(s1n,14)))
    print("   Swrem(da,db): sym %s num %s"%(mp.nstr(ev(s2s),14), mp.nstr(s2n,14)))
