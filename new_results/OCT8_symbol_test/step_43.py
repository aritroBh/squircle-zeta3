import sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
vecfun = vec_of_rational

def S_Li21A(x, y):
    L1 = vecfun(x*y-1); L2 = vecfun(x*y); L3 = vecfun(y-1)
    L4 = vecfun(x*(y-1)/(x-1)); L5 = vecfun(x-1); L6 = vecfun(x)
    return sadd(tensor3(L1,L2,L3), tensor3(L1,L4,L6), tensor3(L3,L5,L6))
def S_Li12_xy(x, y):
    L1 = vecfun(x*y-1); L2 = vecfun(x*(y-1)/(x-1)); L3 = vecfun(y)
    L4 = vecfun(y-1);   L5 = vecfun(x-1);           L6 = vecfun(x*y)
    L7 = vecfun(x/(x-1))
    return sadd(tensor3(L1,L2,L3), tensor3(L4,L5,L3), tensor3(L1,L6,L7), tensor3(L4,L3,L5))
def S_stuffle(a, b):
    # S(Li2(a)) sh S(Li1(b)) = ((1-a)(x)a) sh (1-b)
    T = {}
    av = vecfun(1-a); bv = vecfun(a); cv = vecfun(1-b)
    r = {}
    for (la,ca) in av.items():
        for (lb,cb) in bv.items():
            for (lc,cc) in cv.items():
                for t in [(la,lb,lc),(la,lc,lb),(lc,la,lb)]:
                    r[t] = r.get(t,Fr(0))+ca*cb*cc
    return {k:v for k,v in r.items() if v != 0}

# convB Lam1 = Li12(1/3,-1/3) should equal: S_stuffle(-1/3,1/3) - S_Li21A(-1/3,1/3) - S_Li3(-1/9)
x, y = Fr(1,3), Fr(-1,3)
a, b = y, x   # Li12(b,a): a=-1/3, b=1/3
rhs = sadd(S_stuffle(a,b), smul(S_Li21A(a,b),-1), smul(S_Li3(a*b, vecfun),-1))
lhs = S_Li12_xy(x, y)
diff = sadd(lhs, smul(rhs,-1))
print("S_Li12(1/3,-1/3) vs stuffle-derived:", "MATCH" if not diff else diff)

# also check Lam2 point
x, y = Fr(3,1), Fr(-1,9)
a, b = y, x
rhs = sadd(S_stuffle(a,b), smul(S_Li21A(a,b),-1), smul(S_Li3(a*b, vecfun),-1))
lhs = S_Li12_xy(x, y)
diff = sadd(lhs, smul(rhs,-1))
print("S_Li12(3,-1/9) vs stuffle-derived:", "MATCH" if not diff else diff)
