# step: delta_22(DT_A) and (DT_B) with [2] kept and dropped
import sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
def vsub4(a,b): return tuple(p-q for p,q in zip(a,b))
def vneg4(a): return tuple(-p for p in a)
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
WN={0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
def w2_4(a,b):
    o={}
    for idx,(i,j) in enumerate(B4):
        cc=a[i]*b[j]-a[j]*b[i]
        if cc: o[idx]=o.get(idx,0)+cc
    return o
def spec_delta5_to4(d5, imgi):
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
def spec_sym(s, imgi):
    out={}
    for k,c in s.items():
        nk=tuple(imgi[E.WB5 and v] for v in k) if False else tuple(imgi[v] for v in k)
        out[nk]=out.get(nk,Fr(0))+c
    return {k:v for k,v in out.items() if v}

# specializations
spec_minus = {E.X:vneg4(M), E.U:vneg4(M), E.Y:P, E.V:C, E.W:vsub4(Q,M)}   # y=1-c
spec_plus  = {E.X:vneg4(M), E.U:vneg4(M), E.Y:Q, E.V:C, E.W:vsub4(P,M)}   # y=1+c
def drop2(v4): return v4[1:]  # (c,p,q) 3-dim

def d22_4(sym4):
    out={}
    for (a,b,c,d),cf in sym4.items():
        w={}
        A=w2_4(a,b); B=w2_4(c,d)
        for ia,ca in A.items():
            for ib,cb in B.items():
                if ia==ib: continue
                if ia<ib: w[(ia,ib)]=w.get((ia,ib),Fr(0))+ca*cb
                else: w[(ib,ia)]=w.get((ib,ia),Fr(0))-ca*cb
        for k,v in w.items():
            out[k]=out.get(k,Fr(0))+cf*v
            if out[k]==0: del out[k]
    return out

def d22_4_drop2(sym4):
    out={}
    def w2_3(a,b):
        o={}
        B3=[(i,j) for i in range(3) for j in range(i+1,3)]
        for idx,(i,j) in enumerate(B3):
            cc=a[i]*b[j]-a[j]*b[i]
            if cc: o[idx]=o.get(idx,0)+cc
        return o
    for (a,b,c,d),cf in sym4.items():
        a,b,c,d = drop2(a),drop2(b),drop2(c),drop2(d)
        A=w2_3(a,b); B=w2_3(c,d)
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
WN3={0:'cp',1:'cq',2:'pq'}

s22m = spec_sym(S['Li22'], spec_minus)
s22p = spec_sym(S['Li22'], spec_plus)
s31p = spec_sym(S['Li31'], spec_plus)
s13p = spec_sym(S['Li13'], spec_plus)

print("delta Li22(1/2,1-c) kept:", { (WN[i],WN[j]):str(v) for (i,j),v in d22_4(s22m).items()})
print("delta Li22(1/2,1+c) kept:", { (WN[i],WN[j]):str(v) for (i,j),v in d22_4(s22p).items()})
print("delta Li31(1/2,1+c) kept:", { (WN[i],WN[j]):str(v) for (i,j),v in d22_4(s31p).items()})
print("delta Li13(1/2,1+c) kept:", { (WN[i],WN[j]):str(v) for (i,j),v in d22_4(s13p).items()})
print()
print("delta Li22(1/2,1-c) [2]dropped:", { (WN3[i],WN3[j]):str(v) for (i,j),v in d22_4_drop2(s22m).items()})
print("delta Li22(1/2,1+c) [2]dropped:", { (WN3[i],WN3[j]):str(v) for (i,j),v in d22_4_drop2(s22p).items()})
print("delta Li31(1/2,1+c) [2]dropped:", { (WN3[i],WN3[j]):str(v) for (i,j),v in d22_4_drop2(s31p).items()})
