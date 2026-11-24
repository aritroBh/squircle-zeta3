import sys, pickle
sys.path.insert(0,'<project>/research/08_theory/COR8B_depth_two_two_terms')
sys.path.insert(0,'.')
import cor8b_sym as CB
import cor8br_sym as E
from fractions import Fraction as Fr
import sympy as sp

c = sp.symbols('c')
# pair: x = (1+c)^2/(4c), y = 2/(1+c)
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
def vsub4(a,b): return tuple(p-q for p,q in zip(a,b))
def vscale4(k,a): return tuple(k*p for p in a)
# x = (1+c)^2/(4c):  X = 2Q - 2M - C
X4 = vsub4(vsub4(vscale4(2,Q), vscale4(2,M)), C)
U4 = None
# 1-x = (4c - (1+c)^2)/(4c) = (4c - 1 - 2c - c^2)/(4c) = (-(1-c)^2)/(4c): U = 2P - 2M - C
U4 = vsub4(vsub4(vscale4(2,P), vscale4(2,M)), C)
# y = 2/(1+c): Y = M - Q
Y4 = vsub4(M, Q)
# 1-y = (1+c-2)/(1+c) = (c-1)/(1+c) = -(1-c)/(1+c): V = P - Q
V4 = vsub4(P, Q)
# xy = (1+c)/(2c): 1-xy = (2c-1-c)/(2c) = (c-1)/(2c) = -(1-c)/(2c): W = P - M - C
W4 = vsub4(vsub4(P, M), C)
img_by_index = {0:X4, 1:Y4, 2:U4, 3:V4, 4:W4}
S = E.build_symbols()
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
def w2_4(a,b):
    o={}
    for idx,(i,j) in enumerate(B4):
        cc=a[i]*b[j]-a[j]*b[i]
        if cc: o[idx]=o.get(idx,0)+cc
    return o
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

data = pickle.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/pairs_delta_all_v2.pkl','rb'))
stored22, stored31, stored13 = data[0][2], data[0][3], data[0][4]
for name, stored in [('Li22',stored22), ('Li31',stored31), ('Li13',stored13)]:
    mine = spec_delta5(E.delta22_5(S[name]), img_by_index)
    print(name, "mine :", show(mine))
    print(name, "their db:", show(stored))
    print("   match:", mine==stored)
# stated closed forms at this pair
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
def bneg(w): return {k:-v for k,v in w.items()}
Z4 = vadd4(X4,Y4)
cf22 = wedgebiv(bsum(w2_4(W4,Z4), bneg(w2_4(U4,X4))), w2_4(Y4,V4))
cf31 = bsum(wedgebiv(bneg(w2_4(W4,Z4)), w2_4(Y4,V4)), wedgebiv(w2_4(U4,X4), w2_4(V4,X4)))
cf13 = wedgebiv(w2_4(W4,Z4), w2_4(X4,U4))
print("stated cf22:", show(cf22))
print("stated cf31:", show(cf31))
print("stated cf13:", show(cf13))
print("cf22==db:", cf22==stored22, " cf31==db:", cf31==stored31, " cf13==db:", cf13==stored13)
