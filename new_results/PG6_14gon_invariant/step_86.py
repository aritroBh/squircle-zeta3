import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import mpmath as mp
mp.mp.dps=80
g=S.init()
d=pickle.load(open('delta14.pkl','rb'))
sol_c=d['sol_c']; CRs=d['CRs']
NUM={k:mp.mpf(v) for k,v in pickle.load(open('atoms_num.pkl','rb')).items()}
S14=pickle.load(open('S14.pkl','rb'))
# N3 = sum of the 44 weight-3 atom terms' numerical values
def is_w3(nm): return nm.startswith('L3x') or nm.startswith('Ax') or nm.startswith('Lam') or nm.startswith('T3p') or nm.startswith('T3_')
N3=mp.mpf(0)
for mon,c in S14.items():
    syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
    others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others and is_w3(syms[0]):
        N3 += E.fev(c,80)*NUM[syms[0]]
print("N3 (44 weight-3 atoms) =", mp.nstr(N3,40))
# sum c_z Li3(z)
def fev_ft(cf):
    return E.fev(cf,80)
sumL=mp.mpf(0)
for z,zz in enumerate(CRs):
    c=(sol_c[0][z],sol_c[1][z],sol_c[2][z])
    if any(v!=0 for v in c):
        cf=(c[0],c[1],c[2],1)
        sumL += E.fev(cf,80)*mp.polylog(3, E.fev(zz,80))
print("sum c_z Li3(z) =", mp.nstr(sumL,40))
D = N3 - sumL
print("D = N3 - sumL =", mp.nstr(D,40))
# invisible constants: zeta(3), zeta(2)*log(letter)
z3=mp.zeta(3); z2=mp.pi**2/6
la=mp.log(2*mp.cos(2*mp.pi/7)); lu=mp.log(2*mp.cos(mp.pi/7)); l2=mp.log(2)
basis=[z3, z2*l2, z2*la, z2*lu]
rel=mp.pslq([D]+basis, tol=mp.mpf('1e-70'), maxcoeff=10**12)
print("PSLQ(D; z3,z2l2,z2la,z2lu):", rel)
open('progress.log','a').write('PG6 s_C3: N3 vs sum c_z Li3 diff D; PSLQ invis = %s\n'%(rel is not None))
