# Cross-check: run COR8B's engine (for comparison only) at a generic pair, vs my delta, vs their STATED closed forms
import sys
sys.path.insert(0,'<project>/research/08_theory/COR8B_depth_two_two_terms')
sys.path.insert(0,'.')
import importlib
import cor8b_sym as CB   # their engine, used ONLY to cross-check
import cor8br_sym as E
from fractions import Fraction as Fr

# pick generic letter images: x -> some 4-vector, etc. Use 5 independent images is impossible in 4-dim;
# instead compare in a 5-dim setting by using THEIR wedge funcs generalized... simpler:
# specialize both to a pair (x,y) with distinct non-half values via 4-dim letter space (m,c,p,q)
# e.g. x = c/2, y = 1-c. Then X=[c]-[2]=( -1,1,0,0 ) etc.
S = E.build_symbols()
S22 = S['Li22']; S31 = S['Li31']

# my delta at pair via specializing 5->4
def specialize5to4(s, img):
    out={}
    for k,c in s.items():
        nk=tuple(img[v] for v in k)
        out[nk]=out.get(nk,Fr(0))+c
    return {k:v for k,v in out.items() if v}

M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
def vneg4(a): return tuple(-p for p in a)
# x = c/2, y = 1-c:  X=[c]-[2], U=[1-x]=[(2-c)/2]=[(1+ (1-c))/2] -- NOT an S-unit. Choose x=c, y=1-c:
# X=C, U=P (1-c -> wait 1-x=1-c=P), Y=P, V=C, Z=X+Y=C+P (=[c(1-c)]), W=[1-c(1-c)] not in alphabet -> bad
# choose x=c, y=c: U=V=P, W=[1-c^2]=P+Q
img = {E.X:C, E.Y:C, E.U:P, E.V:P, E.W:vadd4(P,Q), E.Z:vadd4(C,C)}
# their engine: specialize_pair with imgs
imgsCB = {'X':C,'U':P,'Y':C,'V':P,'W':vadd4(P,Q)}
d_mine = E.delta22_5(S22)   # 5-dim
# specialize my 5-dim delta to 4-dim: map basis wedges
WB5 = E.WB5
def spec_delta5(d5, img):
    # map each wedge basis element (i,j) to a 4-dim bivector, then wedge again antisym
    out={}
    B4=[(i,j) for i in range(4) for j in range(i+1,4)]
    def w2(a,b):
        o={}
        for idx,(i,j) in enumerate(B4):
            c=a[i]*b[j]-a[j]*b[i]
            if c: o[idx]=o.get(idx,0)+c
        return o
    for (i,j),v in d5.items():
        a=img[WB5[i][0]]; b=img[WB5[i][1]]
        c2=img[WB5[j][0]]; d2=img[WB5[j][1]]
        A=w2(a,b); B=w2(c2,d2)
        for ia,ca in A.items():
            for ib,cb in B.items():
                if ia==ib: continue
                if ia<ib: out[(ia,ib)]=out.get((ia,ib),Fr(0))+v*ca*cb
                else: out[(ib,ia)]=out.get((ib,ia),Fr(0))-v*ca*cb
    return {k:v for k,v in out.items() if v}

mine_4 = spec_delta5(E.delta22_5(S22), img)
theirs_4 = CB.delta22(CB.specialize_pair(CB.S_Li22(), imgsCB))
WN = {0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
def show(d):
    return " + ".join(f"({v})*[{WN[i]}^{WN[j]}]" for (i,j),v in sorted(d.items())) or "0"
print("delta Li22 at (c,c):  mine(5->4):", show(mine_4))
print("                      COR8B engine  :", show(theirs_4))
print("match:", mine_4==theirs_4)
# now their STATED closed form (W^Z-U^X)^(Y^V), computed at this pair:
# build directly in 4-dim
def w2_4(a,b):
    B4=[(i,j) for i in range(4) for j in range(i+1,4)]
    o={}
    for idx,(i,j) in enumerate(B4):
        cc=a[i]*b[j]-a[j]*b[i]
        if cc: o[idx]=o.get(idx,0)+cc
    return o
def wedgebiv(A,B):
    out={}
    for i,ci in A.items():
        for j,cj in B.items():
            if i==j: continue
            if i<j: out[(i,j)]=out.get((i,j),0)+ci*cj
            else: out[(j,i)]=out.get((j,i),0)-ci*cj
    return {k:v for k,v in out.items() if v}
def bsum(*ws):
    out={}
    for w in ws:
        for k,v in w.items(): out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}
X4,U4,Y4,V4,W4,Z4 = C,P,C,P,vadd4(P,Q),vadd4(C,C)
f1 = bsum(w2_4(W4,Z4), {k:-v for k,v in w2_4(U4,X4).items()})
f2 = w2_4(Y4,V4)
cf = wedgebiv(f1,f2)
print("                      COR8B stated closed form:", show(cf))
