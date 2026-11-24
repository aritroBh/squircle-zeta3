import pickle, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import ftfield as F
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
ORDER=[3,4,1,2,5,6,0]
# degree distribution of term-log monomials
from collections import Counter
deg = Counter(sum(k) for k in E3)
print("E3n term-log monomial degrees:", dict(deg), " total monomials:", len(E3))
# polynomial degrees in t (how many nonzero coeffs)
pdeg = Counter()
for k,v in E3.items():
    nz = sum(1 for vv in v if vv[:3]!=(Fr(0),)*3)
    pdeg[nz]+=1
print("nonzero t-poly coeffs count dist:", dict(pdeg))
# letter vectors for the 7 term-logs' constant parts (from letters14 facs)
L = pickle.load(open('letters14.pkl','rb'))['facs']
print("facs keys:", list(L.keys()))
def lvec(name):
    s,e2,ea,eu = L[name]
    return (e2,ea,eu)  # sign dropped
# term-log constant letter vectors in ORDER [3,4,1,2,5,6,0]
# a3->A0, a4->A1(=A0), a1->B1, a2->B2, a5->B5, a6->B6, a0->h
constvec = {'a3':lvec('A0'),'a4':lvec('A0'),'a1':lvec('B1'),'a2':lvec('B2'),'a5':lvec('B5'),'a6':lvec('B6'),'a0':lvec('h')}
print("const letter vecs:", constvec)
pickle.dump(constvec, open('constvec14.pkl','wb'))
open('progress.log','a').write('PG6 s_B1: E3n degrees %s; const letter vecs extracted\n'%dict(deg))
