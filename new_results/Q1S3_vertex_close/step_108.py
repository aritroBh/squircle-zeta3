import pickle
from fractions import Fraction as Fr
sols=pickle.load(open('sunit_li3.pkl','rb'))
def singleton(vec):
    if len(vec)==1:
        k,c=list(vec.items())[0]
        if c==1: return k
    return None
# index solutions by (vec_z singleton, vec_1z singleton)
idx={}
for i,s in enumerate(sols):
    k1=singleton(s['vec_z']); k2=singleton(s['vec_1z'])
    if k1 and k2: idx[(k1,k2)]=i
# z_{v->w}: vec(z)=w, vec(1-z)=v -> idx[(w,v)]
need=[('a','p2'),('p2','a'),('a','q1'),('q1','a'),('a','q2'),('q2','a')]
for v,w in need:
    i=idx.get((w,v))
    print("z_{%s->%s}: "%(v,w), i if i is not None else "MISSING", sols[i]['z'] if i is not None else "", sols[i]['vec_z'] if i is not None else "", sols[i]['vec_1z'] if i is not None else "")
