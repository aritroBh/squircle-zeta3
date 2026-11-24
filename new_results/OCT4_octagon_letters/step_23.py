import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp, mpmath as mp
from itertools import permutations
mp.mp.dps=30
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
cset=set()
for combo in permutations(list(letters.keys()),4):
    if 't*' not in combo and 'conj' not in combo: continue
    v = sp.simplify(cr(*(letters[k] for k in combo)))
    nv = complex(sp.N(v,25))
    if abs(nv.imag)>1e-9: cset.add(v)
# canonical: pick representative per {z, conj, 1/z, 1/conj} orbit with Im>0, |z|<=1
def canon(v):
    nv = complex(sp.N(v,25))
    cands = [v, sp.conjugate(v), sp.simplify(1/v), sp.simplify(1/sp.conjugate(v))]
    best=None
    for c in cands:
        nc = complex(sp.N(c,25))
        if nc.imag>1e-9 and abs(nc)<=1+1e-9:
            return sp.simplify(c), nc
    # fallback: return the one with positive imag
    for c in cands:
        nc=complex(sp.N(c,25))
        if nc.imag>1e-9: return sp.simplify(c), nc
    return v, nv
canonset={}
for v in cset:
    c,nc = canon(v)
    key=(round(nc.real,12),round(nc.imag,12))
    canonset[key]=(c,nc)
print("Distinct complex cross-ratios (canonical, Im>0, |z|<=1):", len(canonset))
for key in sorted(canonset, key=lambda k:(abs(canonset[k][1]), float(mp.arg(canonset[k][1])))):
    c,nc = canonset[key]
    print("   z =", c, "  |z|=",round(abs(nc),6)," arg/pi=",round(float(mp.arg(nc)/mp.pi),4))
log("step: %d canonical complex cross-ratios (Im>0,|z|<=1), all in Q(zeta_8), args = k*pi/8"%len(canonset))
