import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
t0=time.time()
sd = pickle.load(open(base+"DOD1_dodecagon/S_atom_decomp.pkl","rb"))
coefs = sd['coefs']
r3 = sp.sqrt(3)
l2,l3,lu1 = sp.symbols('l2 l3 lu1')
def pq_of_r(e):
    """e in Q(sqrt3) (no logs) -> (p,q)"""
    e2 = sp.radsimp(sp.expand(e))
    p = sp.simplify((e2 + e2.subs(r3,-r3))/2)
    q = sp.simplify((e2 - e2.subs(r3,-r3))/(2*r3))
    return sp.Rational(p), sp.Rational(q)
def coef_parts(e):
    """return dict: () -> (p,q), ('l2',)->(p,q), ... coefficients in Q(sqrt3)."""
    e = sp.expand(e)
    d = {}
    for mon in [sp.Integer(1), l2, l3, lu1]:
        c = e.coeff(mon) if mon!=1 else sp.expand(e).subs({l2:0,l3:0,lu1:0})
        key = '1' if mon==1 else str(mon)
        if c != 0:
            d[key] = pq_of_r(c)
    # sanity: e == sum
    rec = sum((pq[0]+pq[1]*r3)*(sp.Integer(1) if k=='1' else sp.Symbol(k)) for k,pq in d.items())
    assert sp.simplify(e - rec)==0
    return d
out = {}
for k in coefs:
    out[k] = coef_parts(coefs[k])
print("done, t=%.0fs"%(time.time()-t0))
# check Lam/T3a still match
ac = pickle.load(open(base+"DOD2_dodecagon_finish/atom_coefs_pq.pkl","rb"))
bad = 0
for k in ac:
    p,q = ac[k]
    p2,q2 = out[k].get('1',(0,0))
    if sp.Rational(p)!=p2 or sp.Rational(q)!=q2: bad += 1
print("Lam/T3a match:", bad==0)
from collections import Counter
print(Counter((k[0], tuple(sorted(v.keys()))) for k,v in out.items()))
# display depth-1 families
for fam in ('T2','T3p','Ax','L2x','L3x'):
    ks = [k for k in out if k[0]==fam]
    print("\n==", fam)
    for k in sorted(ks, key=lambda z: z[1]):
        print("  ", k[1], {kk: ("(%s)+(%s)r3"%v) for kk,v in out[k].items()})
pickle.dump(out, open('dod5_atom_coefs_all.pkl','wb'))
print("saved, t=%.0fs"%(time.time()-t0))
