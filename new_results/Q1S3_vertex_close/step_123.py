import pickle
from fractions import Fraction as Fr
sols=pickle.load(open('sunit_li3.pkl','rb'))
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def S_li2(i):
    """-(1-z) x z as dict."""
    out={}
    for k1,c1 in sols[i]['vec_1z'].items():
        for k2,c2 in sols[i]['vec_z'].items():
            out[(k1,k2)]=out.get((k1,k2),Fr(0))-c1*c2
    return out
for (V,W) in [('a','p2'),('a','q1'),('a','q2')]:
    tot={(V,W):Fr(1,2),(W,V):Fr(1,2)}   # sym part from 1/2 log V log W
    for i,c in X[(V,W)]:
        c=Fr(c)
        for k,v in S_li2(i).items():
            tot[k]=tot.get(k,Fr(0))+c*v
    # remove zeros
    tot={k:v for k,v in tot.items() if v}
    print((V,W), "S(X) =", tot, " == target", {(V,W):Fr(1)})
