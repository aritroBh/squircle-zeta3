import sys
sys.path.insert(0,'<project>/research/08_theory/COR8B_depth_two_two_terms')
sys.path.insert(0,'.')
import cor8b_sym as CB
import cor8br_sym as E
from fractions import Fraction as Fr

S = E.build_symbols()
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
img_by_index = {0:C, 1:C, 2:P, 3:P, 4:vadd4(P,Q)}  # X->C,Y->C,U->P,V->P,W->P+Q  (x=c,y=c)
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
def w2_4(a,b):
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
def spec_delta5(d5, imgi):
    out={}
    for (i,j),v in d5.items():
        bi,bj = E.WB5[i]; bk,bl = E.WB5[j]
        A=w2_4(imgi[bi],imgi[bj]); B=w2_4(imgi[bk],imgi[bl])
        for ia,ca in A.items():
            for ib,cb in B.items():
                if ia==ib: continue
                if ia<ib: out[(ia,ib)]=out.get((ia,ib),Fr(0))+v*ca*cb
                else: out[(ib,ia)]=out.get((ib,ia),Fr(0))-v*ca*cb
    return {k:v for k,v in out.items() if v}
WN = {0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
def show(d):
    return " + ".join(f"({v})*[{WN[i]}^{WN[j]}]" for (i,j),v in sorted(d.items())) or "0"

mine_4 = spec_delta5(E.delta22_5(S['Li22']), img_by_index)
imgsCB = {'X':C,'U':P,'Y':C,'V':P,'W':vadd4(P,Q)}
theirs_4 = CB.delta22(CB.specialize_pair(CB.S_Li22(), imgsCB))
print("delta Li22 at (c,c): mine :", show(mine_4))
print("                     CB eng:", show(theirs_4))
print("match:", mine_4==theirs_4)
# stated closed form
X4,U4,Y4,V4,W4,Z4 = C,P,C,P,vadd4(P,Q),vadd4(C,C)
f1 = bsum(w2_4(W4,Z4), {k:-v for k,v in w2_4(U4,X4).items()})
cf = wedgebiv(f1, w2_4(Y4,V4))
print("                     stated:", show(cf))
print("stated == engine:", cf==theirs_4)
