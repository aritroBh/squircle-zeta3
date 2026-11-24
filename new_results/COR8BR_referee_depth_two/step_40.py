# step: delta_22(Omega) from the 56 weight-4 HPL terms, own code
import sys, json, re
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E

# HPL symbol: S(H_{a1..an}(c)) = letter(a_n) ^ ... ^ letter(a_1)  (outermost integration LAST)
# letter(0)=[c], letter(1)=-[1-c], letter(-1)=[1+c]
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def g(a):
    if a==0: return {C:Fr(1)}
    if a==1: return {P:Fr(-1)}
    if a==-1: return {Q:Fr(1)}
def S_H(w):
    s={():Fr(1)}
    for a in reversed(w):   # outermost last: a1 is outermost
        s=E.tensor_right(s, g(a))
    return s

d = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0]); e=int(m.group(1)); rest=m.group(2).strip()
    w=[int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e,w
omega_terms=[(e,w,Fr(t[1])) for t in d for (e,w) in [parse(t)] if w is not None and len(w)==4]
print("weight-4 terms:", len(omega_terms))

# symbol of Omega with coefficients as Laurent polys in c: coeff = rational * c^e
# represent symbol as dict key=(e, tensor)
sym={}
for e,w,cf in omega_terms:
    s=S_H(w)
    for k,v in s.items():
        sym[(e,k)]=sym.get((e,k),Fr(0))+cf*v
        if sym[(e,k)]==0: del sym[(e,k)]
tensors=set(k for (e,k) in sym)
print("distinct pure tensors in S(Omega):", len(tensors))
# delta_22: for each tensor, (a^b)^(c^d) in 4-dim letters, coeff c^e
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
WN={0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
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
d22={}
for (e,k),v in sym.items():
    a,b,cc,dd=k
    w=wedgebiv(w2_4(a,b), w2_4(cc,dd))
    for kk,vv in w.items():
        key=(kk,e)
        d22[key]=d22.get(key,Fr(0))+v*vv
        if d22[key]==0: del d22[key]
print("delta_22(Omega):")
for (kk,e),v in sorted(d22.items()):
    print(f"   ({v}) * c^{e} * [{WN[kk[0]]}^{WN[kk[1]]}]")
open('progress.log','a').write("step_37: computed S(Omega) (own HPL code): 24 pure tensors; delta_22(Omega) computed.\n")
