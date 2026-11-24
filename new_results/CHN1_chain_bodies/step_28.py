import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, symengine as S
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
P8=forms[(8,3)]  # zeta(3,5)
print('P8 #monomials:', len(P8))
# pilot point checks
def pe(p,a): return S.peval(p,a)
for a,want in [((1,1,1,1),F(0)),((1,1,1,2),F(58)),((1,2,2,1),F(-1464,5)),((3,1,2,5),F(36696)),((2,1,1,1),F(58)),((5,2,1,3),F(36696))]:
    got=pe(P8,a); print(a, got, 'PASS' if got==want else 'FAIL want %s'%want)
# divisibility by a1*a2*a3*a4
div = all(all(x>=1 for x in e) for e in P8)
print('P8 divisible by a1a2a3a4:', div)
Q4={tuple(x-1 for x in e):c for e,c in P8.items()}
print('Q4 #monomials:', len(Q4))
# print Q4 nicely sorted
def pstr(p):
    vs='xyzw'
    terms=[]
    for e in sorted(p, key=lambda e:(-sum(e),e)):
        c=p[e]
        m=''.join(vs[i]+('' if e[i]==1 else '^%d'%e[i]) for i in range(4) if e[i]>0)
        terms.append('%+s%s'%(c, '*'+m if m else ''))
    return ' '.join(terms)
print('Q4 =', pstr(Q4))
print('Q4(1,1,1,1) =', pe(Q4,(1,1,1,1)))
print('t=%.1fs'%(time.time()-t0))
