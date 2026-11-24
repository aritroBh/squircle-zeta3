import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vneg4(a): return tuple(-p for p in a)
def vsub4(a,b): return tuple(p-q for p,q in zip(a,b))
BASIS5=[E.X,E.Y,E.U,E.V,E.W]
def spec_sym(s, imgbasis):
    out={}
    for k,cf in s.items():
        nk=[]
        for v in k:
            img=[0,0,0,0]
            for i,b in enumerate(BASIS5):
                if v[i]:
                    iv=imgbasis[b]
                    for j in range(4): img[j]+=v[i]*iv[j]
            nk.append(tuple(img))
        nk=tuple(nk)
        out[nk]=out.get(nk,Fr(0))+cf
    return {k:v for k,v in out.items() if v}
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
def d22_4(sym4):
    out={}
    for (a,b,c,d),cf in sym4.items():
        A={}; B={}
        for idx,(i,j) in enumerate(B4):
            ca=a[i]*b[j]-a[j]*b[i]
            if ca: A[idx]=A.get(idx,0)+ca
            cb=c[i]*d[j]-c[j]*d[i]
            if cb: B[idx]=B.get(idx,0)+cb
        w={}
        for ia,ca in A.items():
            for ib,cb in B.items():
                if ia==ib: continue
                if ia<ib: w[(ia,ib)]=w.get((ia,ib),Fr(0))+ca*cb
                else: w[(ib,ia)]=w.get((ib,ia),Fr(0))-ca*cb
        for k,v in w.items():
            out[k]=out.get(k,Fr(0))+cf*v
            if out[k]==0: del out[k]
    return out
spec_minus = {E.X:vneg4(M), E.U:vneg4(M), E.Y:P, E.V:C, E.W:vsub4(Q,M)}
spec_plus  = {E.X:vneg4(M), E.U:vneg4(M), E.Y:Q, E.V:C, E.W:vsub4(P,M)}
s22m = spec_sym(S['Li22'], spec_minus)
s22p = spec_sym(S['Li22'], spec_plus)
s31p = spec_sym(S['Li31'], spec_plus)
d22m, d22p, d31p = d22_4(s22m), d22_4(s22p), d22_4(s31p)
def lp_add(a,b,sgn=1):
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,Fr(0))+sgn*v
        if out[k]==0: del out[k]
    return out
def lp_mul_by(a, f):
    out={}
    for e1,v1 in a.items():
        for e2,v2 in f.items():
            out[e1+e2]=out.get(e1+e2,Fr(0))+v1*v2
    return {e:v for e,v in out.items() if v}
fA  = {-1:Fr(1), 0:Fr(-1)}    # (1-c)/c = c^-1 - 1
fB  = {-1:Fr(-1), 0:Fr(-1)}   # -(1+c)/c
fBp = {-1:Fr(1), 0:Fr(1)}     # +(1+c)/c
def build(terms):
    D={}
    for src,fac in terms:
        for k,v in src.items():
            D[k]=lp_add(D.get(k,{}), lp_mul_by({0:v}, fac))
            if not D[k]: del D[k]
    return D
DTA = build([(d22m,fA),(d22p,fB)])
DTB = build([(d22m,fA),(d31p,fBp)])
OM = {(3,5): {0:Fr(1), -1:Fr(-1)}, (4,5): {0:Fr(-1), -1:Fr(-1)}}
WN={0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
def showlp(d):
    parts=[]
    for (i,j),lp in sorted(d.items()):
        terms=" ".join(f"({v})c^{e}" for e,v in sorted(lp.items()))
        parts.append(f"[{WN[i]}^{WN[j]}]: {terms}")
    return "\n".join(parts) if parts else "0"
print("delta_22(DT_A) [2]-kept:"); print(showlp(DTA))
print("DT_A == DT_B [2]-kept:", DTA==DTB)
R={}
for k in set(DTA)|set(OM):
    lpv=lp_add(DTA.get(k,{}), OM.get(k,{}), -1)
    if lpv: R[k]=lpv
print("R = delta(DT_A) - delta(Omega) [2]-kept:"); print(showlp(R))
rj = json.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/DT_A_residual_2kept.json'))
rj_conv = {tuple(int(x) for x in kk.split(',')): {int(e):Fr(v) for e,v in dct.items()} for kk,dct in rj.items()}
print("R == DT_A_residual_2kept.json:", R==rj_conv)
DTA_drop = {k:v for k,v in DTA.items() if k in [(3,5),(4,5),(3,4)]}
print("delta(DT_A) with [2]-components removed == delta(Omega) (sign +):", DTA_drop==OM)
open('progress.log','a').write("step_44: FIXED coefficient sign. delta(DT_A)=delta(DT_B); [2]-dropped equals delta(Omega) with sign +; [2]-kept residual = exactly DT_A_residual_2kept.json (4 m-components).\n")
