import sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g=S.init()
def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items(): r[k]=r.get(k,Fr(0))+c
    return {k:c for k,c in r.items() if c}
def tensor(u,v,w,coeff=Fr(1)):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            for c,cc in w.items(): out[(a,b,c)]=out.get((a,b,c),Fr(0))+coeff*ca*cb*cc
    return {k:x for k,x in out.items() if x}
def sadd(*Ss):
    r={}
    for s in Ss:
        for k,c in s.items(): r[k]=r.get(k,Fr(0))+c
    return {k:c for k,c in r.items() if c}
def sscale(s,q):
    q=Fr(q); return {k:c*q for k,c in s.items() if c*q}
# S_Li21(x,y) my convention (from sym_Lam core)
def S_Li21(x,y):
    L1=S.vec(E.fsub(E.fmul(x,y),E.ONE)); L2=S.vec(E.fmul(x,y)); L3=S.vec(E.fsub(y,E.ONE))
    L4=S.vec(E.fdiv(E.fmul(x,E.fsub(y,E.ONE)),E.fsub(x,E.ONE))); L5=S.vec(E.fsub(x,E.ONE)); L6=S.vec(x)
    Sdod=sadd(tensor(L1,L2,L3),tensor(L1,L4,L6),tensor(L3,L5,L6))
    out={}
    for (a,b,c),cf in Sdod.items():
        out[(c,b,a)]=out.get((c,b,a),Fr(0))-cf
        if out[(c,b,a)]==0: del out[(c,b,a)]
    return out
# S_Li12(x,y) via G-word I(0;1/(xy),0,1/y;1) using S3 (first letter nonzero)
def S_Li12(x,y):
    a=E.finv(E.fmul(x,y)); b=E.ZERO; c=E.finv(y)
    return S.S3(a,b,c)
def S_Li3(z): return S.li3_symbol(z)
def S_Li2(x):  # x (x) (1-x), weight-2 analog of li3 convention
    zv=S.vec(x); omz=S.vec(E.fsub(E.ONE,x))
    return tensor(zv,omz,{} if False else {0:Fr(0)})  # placeholder; handle separately
def shuffle_Li2_Li1(x,y, sgn=1):
    # S(Li2(x)) [w2] III S(Li1(y)) [w1]: shuffle of 2-tensor and 1-tensor -> 3-tensor
    zv=S.vec(x); omz=S.vec(E.fsub(E.ONE,x)); yv=S.vec(E.fsub(E.ONE,y))  # S(Li1(y)) = (1-y) * sgn?
    w2=tensor(zv,omz,{0:Fr(1)})  # hack: build (a,b) pairs
    # build 2-tensor x (x) (1-x)
    T2={}
    for a,ca in zv.items():
        for b,cb in omz.items(): T2[(a,b)]=T2.get((a,b),Fr(0))+ca*cb
    out={}
    for (a,b),cab in T2.items():
        for c,cc in yv.items():
            # shuffle: c can go in 3 positions
            for key in [(c,a,b),(a,c,b),(a,b,c)]:
                out[key]=out.get(key,Fr(0))+sgn*cab*cc
    return {k:v for k,v in out.items() if v}
x=E.X[1]; y=E.fdiv(E.X[2],E.X[1])
xy=E.fmul(x,y)
LHS1 = S_Li21(x,y)
LHS2 = S_Li12(y,x)
LHS3 = S_Li3(xy)
LHS = sadd(LHS1,LHS2,LHS3)
for sgn in [1,-1]:
    RHS = shuffle_Li2_Li1(x,y,sgn)
    R=dict(LHS)
    for k,c in RHS.items():
        R[k]=R.get(k,Fr(0))-c
        if R[k]==0: del R[k]
    print(f"stuffle residual (Li1 sign {sgn}): {len(R)} terms", R if len(R)<12 else "")
