import sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
def vsub4(a,b): return tuple(p-q for p,q in zip(a,b))
def vneg4(a): return tuple(-p for p in a)
BASIS5=[E.X,E.Y,E.U,E.V,E.W]
def spec_sym(s, imgbasis):
    # imgbasis: dict basis-vec -> 4-vector; letters may be linear combos
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
B3=[(i,j) for i in range(3) for j in range(i+1,3)]
WN={0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
WN3={0:'cp',1:'cq',2:'pq'}
def wedgebiv_from(A,B, WB):
    w={}
    for ia,ca in A.items():
        for ib,cb in B.items():
            if ia==ib: continue
            if ia<ib: w[(ia,ib)]=w.get((ia,ib),Fr(0))+ca*cb
            else: w[(ib,ia)]=w.get((ib,ia),Fr(0))-ca*cb
    return {k:v for k,v in w.items() if v}
def d22_4(sym4):
    out={}
    for (a,b,c,d),cf in sym4.items():
        A={}; B={}
        for idx,(i,j) in enumerate(B4):
            ca=a[i]*b[j]-a[j]*b[i]
            if ca: A[idx]=A.get(idx,0)+ca
            cb=c[i]*d[j]-c[j]*d[i]
            if cb: B[idx]=B.get(idx,0)+cb
        for k,v in wedgebiv_from(A,B,B4).items():
            out[k]=out.get(k,Fr(0))+cf*v
            if out[k]==0: del out[k]
    return out
def d22_drop2(sym4):
    out={}
    for (a,b,c,d),cf in sym4.items():
        a,b,c,d=a[1:],b[1:],c[1:],d[1:]
        A={}; B={}
        for idx,(i,j) in enumerate(B3):
            ca=a[i]*b[j]-a[j]*b[i]
            if ca: A[idx]=A.get(idx,0)+ca
            cb=c[i]*d[j]-c[j]*d[i]
            if cb: B[idx]=B.get(idx,0)+cb
        for k,v in wedgebiv_from(A,B,B3).items():
            out[k]=out.get(k,Fr(0))+cf*v
            if out[k]==0: del out[k]
    return out
spec_minus = {E.X:vneg4(M), E.U:vneg4(M), E.Y:P, E.V:C, E.W:vsub4(Q,M)}
spec_plus  = {E.X:vneg4(M), E.U:vneg4(M), E.Y:Q, E.V:C, E.W:vsub4(P,M)}
s22m = spec_sym(S['Li22'], spec_minus)
s22p = spec_sym(S['Li22'], spec_plus)
s31p = spec_sym(S['Li31'], spec_plus)
s13p = spec_sym(S['Li13'], spec_plus)
def show4(d): return " + ".join(f"({v})*[{WN[i]}^{WN[j]}]" for (i,j),v in sorted(d.items())) or "0"
def show3(d): return " + ".join(f"({v})*[{WN3[i]}^{WN3[j]}]" for (i,j),v in sorted(d.items())) or "0"
print("d Li22(1/2,1-c) kept:", show4(d22_4(s22m)))
print("d Li22(1/2,1+c) kept:", show4(d22_4(s22p)))
print("d Li31(1/2,1+c) kept:", show4(d22_4(s31p)))
print("d Li13(1/2,1+c) kept:", show4(d22_4(s13p)))
print()
print("d Li22(1/2,1-c) dropped:", show3(d22_drop2(s22m)))
print("d Li22(1/2,1+c) dropped:", show3(d22_drop2(s22p)))
print("d Li31(1/2,1+c) dropped:", show3(d22_drop2(s31p)))
import pickle
pickle.dump({'s22m':s22m,'s22p':s22p,'s31p':s31p,'s13p':s13p}, open('my_dt_syms.pkl','wb'))
