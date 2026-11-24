# step_19: delta(N3) for the regular octagon, decagon, 12-gon
from fractions import Fraction as Fr
import pickle

def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def vmul(v,a):
    a=Fr(a); return {k:c*a for k,c in v.items() if c*a}
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            if a<b: key,s=(a,b),1
            else: key,s=(b,a),-1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3v(terms):
    """terms: (p, q, vz, vomz): coeff = p + q*sqrt(D); delta summed as pair (P, Q)."""
    P, Q = {}, {}
    for p,q,vz,vomz in terms:
        w = wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                if p:
                    P[key]=P.get(key,Fr(0))+p*cab*c3
                    if P[key]==0: del P[key]
                if q:
                    Q[key]=Q.get(key,Fr(0))+q*cab*c3
                    if Q[key]==0: del Q[key]
    return P,Q

T = pickle.load(open("field_tables.pkl","rb"))
T2,T5,T3 = T["T2"],T["T5"],T["T3"]

# ---- regular octagon ----
oct_terms = []
for nm, (p,q) in [("d",(Fr(3,8),Fr(-1,4))), ("k",(Fr(-3,8),Fr(1,4))),
                  ("mr2h",(Fr(-9,2),Fr(3))), ("r2h",(Fr(-155,32),Fr(103,32))),
                  ("x1",(Fr(33,32),Fr(-21,32))), ("x2",(Fr(-79,32),Fr(59,32)))]:
    vz = T2[nm][1]; vomz = T2["1-"+nm][1]
    oct_terms.append((p,q,vz,vomz))
P,Q = delta_li3v(oct_terms)
print("=== regular octagon: delta(N3) = (P part, sqrt2 part)")
print(" P:", P if P else 0)
print(" Q:", Q if Q else 0)
print(" delta == 0 :", not P and not Q)

# ---- decagon ----
dec_terms = []
for nm,(p,q) in [("1/phi",(Fr(-75,32),Fr(81,64))), ("-1/phi",(Fr(-35,32),Fr(53,64)))]:
    vz = T5[nm][1]
    omz = "1-1/phi" if nm=="1/phi" else "1+1/phi"
    vomz = T5[omz][1]
    dec_terms.append((p,q,vz,vomz))
P,Q = delta_li3v(dec_terms)
print("=== decagon: delta(N3)")
print(" P:", P if P else 0, " Q:", Q if Q else 0, " delta==0:", not P and not Q)

# ---- 12-gon ----
dod_terms = []
coeffs = [("z1",Fr(-73,24),Fr(37,8)), ("z2",Fr(11,2),Fr(-39,4)), ("z3",Fr(97,24),Fr(-63,8)),
          ("z4",Fr(19,192),Fr(-13,64)), ("z5",Fr(-1,8),Fr(9,8)), ("z6",Fr(2,3),Fr(-3,4)),
          ("z7",Fr(7,4),Fr(-13,4)), ("z8",Fr(-17,96),Fr(13,16))]
for nm,p,q in coeffs:
    vz = T3[nm][1]; vomz = T3["1-"+nm][1]
    dod_terms.append((p,q,vz,vomz))
P,Q = delta_li3v(dod_terms)
print("=== 12-gon: delta(N3)")
print(" P:", P if P else 0)
print(" Q:", Q if Q else 0)
print(" delta == 0 :", not P and not Q)
