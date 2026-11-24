import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
from itertools import permutations
sqrt2 = sp.sqrt(2); I = sp.I; oo = sp.oo
tstar = (1+I*(1+sqrt2))/2
letters = {'0':sp.Integer(0),'1':sp.Integer(1),'inf':oo,
           '-(1+r2)':-(1+sqrt2),'2+r2':2+sqrt2,'t*':tstar,'conj':sp.conjugate(tstar)}
def term(x,y):
    if x is oo and y is oo: return sp.Integer(1)
    if x is oo: return sp.Integer(1)
    if y is oo: return sp.Integer(-1)
    return x-y
def cr(a,b,c,d):
    def ratio(p,q,r,ss):
        num = term(p,q); den = term(r,ss)
        if num==1 and den==1: return sp.Integer(1)
        return num/den
    return sp.simplify(ratio(a,c,a,d)/ratio(b,c,b,d))
allkeys = list(letters.keys())
cset={}
for combo in permutations(allkeys,4):
    if 't*' not in combo and 'conj' not in combo: continue
    v = sp.simplify(cr(*(letters[k] for k in combo)))
    cset[sp.nsimplify(v)] = None
print("Cross-ratios involving complex letters:", len(cset))
reals=[]; complexes=[]
for v in cset:
    nv = complex(sp.N(v,20))
    (complexes if abs(nv.imag)>1e-9 else reals).append((v,nv))
print("real-valued:", len(reals), " genuinely complex:", len(complexes))
print("\n--- real ones ---")
for v,nv in sorted(reals,key=lambda z:z[1].real):
    print("   ", v, "=", round(nv.real,6))
print("\n--- complex ones (value, |z|, arg/pi) ---")
import mpmath as mp; mp.mp.dps=30
for v,nv in sorted(complexes,key=lambda z:(abs(z[1]),z[1].imag))[:40]:
    print("   ", v, "  |z|=",round(abs(nv),5), " arg/pi=", round(float(mp.arg(nv)/mp.pi),5))
log("step: complex-letter cross-ratios: %d real, %d complex"%(len(reals),len(complexes)))
