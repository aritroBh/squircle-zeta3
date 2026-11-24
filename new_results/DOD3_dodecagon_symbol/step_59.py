import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *
def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d
r3 = mp.sqrt(3)
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+r3); z2,z3c = mp.pi**2/6, mp.zeta(3)
lp2, lp3, le = (l2+lu)/2, l3/2, lu
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
bidx_zs = [(2,-1,1),(-1,1,2),(0,1,3),(7,-4,1),(2,-1,2),(-3,2,3),(3,-1,6),(1,0,3)]
bidx_zs = [P(*z) for z in bidx_zs]
vals = [mp.re(mp.polylog(3, mp.mpc(mpf_of(z0)))) for z0 in [P(1,0,4)]+bidx_zs]
names = ['Li3(1/4)'] + ['Li3(%s)'%str(z) for z in bidx_zs]
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        vals.append(lp2**a*lp3**b*le**c); names.append("lp2^%d lp3^%d le^%d"%(a,b,c))
for w in b2zs:
    for p in ('2','3','e'):
        vals.append(mp.polylog(2, mpf_of(w))*{'2':lp2,'3':lp3,'e':le}[p])
        names.append("Li2(%s)*l%s"%(str(w),p))
vals += [z3c, z2*l2, z2*l3, z2*lu, L3c]
names += ['z3','z2l2','z2l3','z2lu','L3chi']
rel = mp.pslq(vals, tol=mp.mpf('1e-550'), maxcoeff=10**18, maxsteps=300000)
print("PSLQ result:", rel)
if rel:
    resid = abs(sum(mp.mpf(r)*v for r,v in zip(rel, vals)))
    print("resid:", mp.nstr(resid,3), " maxc:", max(abs(r) for r in rel))
    for n,r in zip(names,rel):
        if r: print("   %-24s %s"%(n,r))
