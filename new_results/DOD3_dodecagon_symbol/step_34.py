import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 200
import sympy as sp
from dod3_R import *
base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon'); d2 = os.path.join(base,'DOD2_dodecagon_finish')
atoms620 = pickle.load(open(os.path.join(d1,'atoms620.pkl'),'rb'))
g = pickle.load(open(os.path.join(d2,'gdswaps620.pkl'),'rb'))

vals200 = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals200)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,300))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,300)))))
    return mp.mpf(sp.N(e, 200))

for (i,j,k) in [(6,1,5),(1,2,5),(1,2,6),(2,5,6)]:
    xi = Xs[i]
    R = build_R(i,j,k)
    rv = ev(R)
    rq = atoms620['T3a'][(i,j,k)] - (g['SwG'][(i,j,k)]-g['SwD'][(i,j,k)])/atoms620['X'][i]
    print((i,j,k), " R_sym =", mp.nstr(rv,16))
    print("          R_quad=", mp.nstr(rq,16), "  diff =", mp.nstr(abs(rv-rq),3))
