import sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g=S.init()
def tensor(u,v,w,coeff=Fr(1)):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            for c,cc in (w.items() if w else [(None,Fr(1))]):
                key=(a,b) if w is None else (a,b,c)
                out[key]=out.get(key,Fr(0))+coeff*ca*cb*cc
    return {k:x for k,x in out.items() if x}
def sadd(*Ss):
    r={}
    for s in Ss:
        for k,c in s.items(): r[k]=r.get(k,Fr(0))+c
    return {k:c for k,c in r.items() if c}
def sscale(s,q):
    q=Fr(q); return {k:c*q for k,c in s.items() if c*q}
def S_Li21(x,y):
    L1=S.vec(E.fsub(E.fmul(x,y),E.ONE)); L2=S.vec(E.fmul(x,y)); L3=S.vec(E.fsub(y,E.ONE))
    L4=S.vec(E.fdiv(E.fmul(x,E.fsub(y,E.ONE)),E.fsub(x,E.ONE))); L5=S.vec(E.fsub(x,E.ONE)); L6=S.vec(x)
    Sdod=sadd(tensor(L1,L2,L3),tensor(L1,L4,L6),tensor(L3,L5,L6))
    out={}
    for (a,b,c),cf in Sdod.items():
        out[(c,b,a)]=out.get((c,b,a),Fr(0))-cf
        if out[(c,b,a)]==0: del out[(c,b,a)]
    return out
def S_Li3(z): return S.li3_symbol(z)
# S(Li2(x)) = x (x) (1-x)  [2-tensor]
def S_Li2(x):
    zv=S.vec(x); omz=S.vec(E.fsub(E.ONE,x)); T={}
    for a,ca in zv.items():
        for b,cb in omz.items(): T[(a,b)]=T.get((a,b),Fr(0))+ca*cb
    return T
# S(Li1(y)) = -(1-y)  [1-tensor]  (Li1 = -log(1-y))
def S_Li1(y):
    yv=S.vec(E.fsub(E.ONE,y)); return { (k,):-c for k,c in yv.items() if c}
def shuffle21(T2,T1):
    # shuffle a 2-tensor and a 1-tensor -> 3-tensor (insert the 1-letter in 3 positions)
    out={}
    for (a,b),cab in T2.items():
        for (c,),cc in T1.items():
            for key in [(c,a,b),(a,c,b),(a,b,c)]:
                out[key]=out.get(key,Fr(0))+cab*cc
    return {k:v for k,v in out.items() if v}
# S(Li_{1,2}(y,x)) via G-word I(0;1/(xy),0,1/x;1)
def S_Li12_yx(x,y):
    a=E.finv(E.fmul(x,y)); c=E.finv(x)
    return S.S3(a,E.ZERO,c)
x=E.X[1]; y=E.fdiv(E.X[2],E.X[1]); xy=E.fmul(x,y)
LHS=sadd(S_Li21(x,y), S_Li12_yx(x,y), S_Li3(xy))
RHS=shuffle21(S_Li2(x), S_Li1(y))
R=dict(LHS)
for k,c in RHS.items():
    R[k]=R.get(k,Fr(0))-c
    if R[k]==0: del R[k]
print("stuffle residual:", len(R), R if len(R)<20 else "(many)")
