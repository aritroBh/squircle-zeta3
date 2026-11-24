# chkD_p3.py -- the p=3 instance of E3 / adjacent identities / diagonal polynomial, exact, with finiteness audit.
import sys, pickle, os
from fractions import Fraction as F
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import prove_E3b56 as E3
import prove_adj56 as AD
from math import comb
p=E3.p
R='<project>/research/08_theory/'
H=pickle.load(open(R+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def finite(x): return not (x.has(sp.zoo) or x.has(sp.nan) or x.has(sp.oo))
# 1. Delta_pp(delta) at p=3: exceptional offset of the diagonal (3,3) at k=6, per depth
print('Delta_33(delta)/rho_33(0)-type offsets (exact, p=3):')
for delta in range(0,7):
    D=E3.Delta(delta)
    v=sp.nsimplify(D.subs(p,3))
    print(f'  delta={delta}: Delta(p=3) = {v}', ' finite' if finite(v) else ' POLE')
# 2. E3 feed at p=3,4,5,6,7 for delta=0..6
print('E3 feed into child (p,p-1) at p=3..7 (should be 0 for delta<=6):')
for delta in range(0,7):
    f=E3.feed(delta)
    vals=[sp.simplify(f.subs(p,pp)) for pp in (3,4,5,6,7)]
    print(f'  delta={delta}: {vals}')
# 3. diagonal polynomial at p=3 vs exact data (3,3): u(delta)(k)/rho_33(0) == sum_s C(k-7,s) ehat(delta,s)
def to_frac_pair(expr):
    expr=sp.expand(sp.nsimplify(expr)); re_,im_=expr.as_real_imag()
    return F(int(re_.p),int(re_.q)),F(int(im_.p),int(im_.q))
rho0=H[(3,3)][0][6]   # (p0,p1,p2,p3)
sys.path.insert(0, R+'WTS_weight_k_structure/theory')
from fastred import QSI
r0=QSI(*rho0)
bad=0;n=0
EH={}
for delta in range(0,7):
    for s in range(0,delta//2+1):
        EH[(delta,s)]=to_frac_pair(E3.ehat(delta,s).subs(p,3))
for delta in range(0,7):
    for k in range(7,53):
        x=F(0);y=F(0)
        for s in range(0,delta//2+1):
            ex,ey=EH[(delta,s)]; c=comb(k-7,s); x+=c*ex; y+=c*ey
        pred=r0*QSI(x,0,y,0)
        n+=1
        if (pred.p0,pred.p1,pred.p2,pred.p3)!=tuple(H[(3,3)][delta][k]): bad+=1; print('diag poly mismatch',delta,k)
print(f'diagonal (3,3) stable polynomial (ehat at p=3) vs exact towers: {n} cells (delta 0..6, k=7..52), {bad} mismatches')
# 4. adjacent identities at p=3 and 4 exactly, with finiteness audit of every term
import itertools
allok=True
for (delta,s) in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    terms=AD.feed_terms(delta,s)
    lhs=sp.expand(AD.eR(delta,s,AD.d,p)+AD.eR(delta,s+1,AD.d,p))
    for pp in (3,4,5):
        tv=[sp.simplify(t.subs(p,pp)) for t in terms]
        lv=sp.simplify(lhs.subs(p,pp))
        fin=all(finite(t) for t in tv) and finite(lv)
        res=sp.simplify(lv-sum(tv))
        ok=(res==0) and fin
        allok&=ok
        if pp==3 or not ok: print(f'  adjacent ({delta},{s}) at p={pp}: residual {res}, all terms finite: {fin}')
print('adjacent identities at p=3,4,5 exact, all finite:',allok)
