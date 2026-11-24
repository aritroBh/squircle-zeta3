import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import mpmath as mp
mp.mp.dps=70
g=S.init()
S14=pickle.load(open('S14.pkl','rb'))
mons14=pickle.load(open('mons14.pkl','rb'))
# Recompute z3 coefficient tagged by source: local (wc=0), single (1 fam), pair (2), triple (3)
IDX=[1,2,5,6]
INTC={}
def INTc(n,a,b,wc):
    k=(n,a,b,wc)
    if k not in INTC: INTC[k]=E.INT(n,a,b,wc)
    return INTC[k]
z3_by_source={0:E.ZERO,1:E.ZERO,2:E.ZERO,3:E.ZERO}
for (n,ab,wc,fk),coef in mons14.items():
    nfam=sum(1 for c in wc if c>0)
    J=INTc(n,ab[0],ab[1],wc)
    # letter monomial
    lm={}
    if fk[0]: lm['l2']=fk[0]
    if fk[1]: lm['la']=fk[1]
    if fk[2]: lm['lu']=fk[2]
    lm=E.mkey(lm)
    for mon,c in J.items():
        d=dict(mon)
        if d.get('z3')==1 and all(s in ('z3','l2','la','lu') for s in d):
            # z3 coefficient (with letter monomial)
            z3_by_source[nfam]=E.fadd(z3_by_source[nfam], E.fmul(c,coef))
LL=mp.sqrt(E.fev(E._g['L'],70))
C7 = 48*7*mp.sin(mp.pi/14)**4*mp.cos(mp.pi/14)*(3*mp.sin(mp.pi/14)**2-1)
print("z3 coeff by #far-families:")
for k in range(4):
    print(f"  nfam={k}: {E.fev(z3_by_source[k],30)}")
loc = E.fev(z3_by_source[0],70)
print("84L * (local z3, nfam=0) =", mp.nstr(84*LL*loc,30))
print("C_7                          =", mp.nstr(C7,30))
print("diff:", mp.nstr(abs(84*LL*loc-C7),5))
# also nfam=0+1 (singles)
locs = E.fev(E.fadd(z3_by_source[0],z3_by_source[1]),70)
print("84L * (local+single z3) =", mp.nstr(84*LL*locs,30))
open('progress.log','a').write('PG6 s_C2: C7 check: 84L*local z3 = %s vs C7 %s\n'%(mp.nstr(84*LL*loc,20),mp.nstr(C7,20)))
