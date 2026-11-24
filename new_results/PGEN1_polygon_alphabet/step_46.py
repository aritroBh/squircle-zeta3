# step_22: (i) delta^2=0 (wedge3) for all polygons; (ii) 12-gon delta invariance under duplication z4=z1^2
from fractions import Fraction as Fr
import pickle

def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
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
    P,Q={},{}
    for p,q,vz,vomz in terms:
        w=wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                if p: P[key]=P.get(key,Fr(0))+p*cab*c3
                if q: Q[key]=Q.get(key,Fr(0))+q*cab*c3
    return ({k:v for k,v in P.items() if v},{k:v for k,v in Q.items() if v})
def wedge3_image(D):
    """map ((a,b),c) -> a^b^c ; return dict (a,b,c)->coef"""
    out={}
    for ((a,b),c),coef in D.items():
        letters = sorted(set([a,b,c]))
        if len(letters)<3: continue
        # sign of permutation from (a,b,c) to sorted
        perm = [a,b,c]; srt = sorted(perm)
        inv = 0
        for i in range(3):
            for j in range(i+1,3):
                if perm[i]>perm[j]: inv+=1
        sgn = -1 if inv%2 else 1
        key = tuple(srt)
        out[key]=out.get(key,Fr(0))+sgn*coef
    return {k:v for k,v in out.items() if v}

T = pickle.load(open("field_tables.pkl","rb"))
T2,T5,T3 = T["T2"],T["T5"],T["T3"]

# octagon
oct_terms = [(p,q,T2[nm][1],T2["1-"+nm][1]) for nm,(p,q) in
    [("d",(Fr(3,8),Fr(-1,4))),("k",(Fr(-3,8),Fr(1,4))),("mr2h",(Fr(-9,2),Fr(3))),
     ("r2h",(Fr(-155,32),Fr(103,32))),("x1",(Fr(33,32),Fr(-21,32))),("x2",(Fr(-79,32),Fr(59,32)))]]
P,Q = delta_li3v(oct_terms)
print("octagon wedge3(delta):", wedge3_image(P), wedge3_image(Q), " (empty => in ker, i.e. in Lie^3)")

# 12-gon standard
coeffs = [("z1",Fr(-73,24),Fr(37,8)),("z2",Fr(11,2),Fr(-39,4)),("z3",Fr(97,24),Fr(-63,8)),
          ("z4",Fr(19,192),Fr(-13,64)),("z5",Fr(-1,8),Fr(9,8)),("z6",Fr(2,3),Fr(-3,4)),
          ("z7",Fr(7,4),Fr(-13,4)),("z8",Fr(-17,96),Fr(13,16))]
dod_terms = [(p,q,T3[nm][1],T3["1-"+nm][1]) for nm,p,q in coeffs]
P,Q = delta_li3v(dod_terms)
print("12-gon wedge3(delta):", wedge3_image(P), wedge3_image(Q))

# (ii) duplicate z4 -> 4 z1 + 4 (-z1): vec(-z1)=vec(z1); vec(1+z1) = vec(3-sqrt3) = g2*g3/u
vz1 = T3["z1"][1]
v1pz1 = {"g2":Fr(1),"g3":Fr(1),"u":Fr(-1)}
import sympy as sp
assert sp.simplify((1+(2-sp.sqrt(3))) - (1+sp.sqrt(3))*sp.sqrt(3)/(2+sp.sqrt(3)))==0
dod_terms2 = []
for nm,p,q in coeffs:
    if nm=="z4":
        dod_terms2.append((4*p,4*q,T3["z1"][1],T3["1-z1"][1]))
        dod_terms2.append((4*p,4*q,T3["z1"][1],v1pz1))  # {-z1}: vec same, 1-(-z1)=1+z1
    else:
        dod_terms2.append((p,q,T3[nm][1],T3["1-"+nm][1]))
P2,Q2 = delta_li3v(dod_terms2)
def diff(A,B): return vadd(A, {k:-c for k,c in B.items()})
print("12-gon delta: 8-pt vs duplicated(z4->4z1+4(-z1)) presentation identical:",
      diff(P,P2)=={} and diff(Q,Q2)=={})
# decagon wedge3 trivially 0; print delta again
dec_terms = []
for nm,(p,q),omz in [("1/phi",(Fr(-75,32),Fr(81,64)),"1-1/phi"), ("-1/phi",(Fr(-35,32),Fr(53,64)),"1+1/phi")]:
    dec_terms.append((p,q,T5[nm][1],T5[omz][1]))
P3,Q3 = delta_li3v(dec_terms)
print("decagon delta:", P3, Q3)
