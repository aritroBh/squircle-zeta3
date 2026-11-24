# step_20: debug octagon delta term by term
from fractions import Fraction as Fr
import pickle
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            if a<b: key,s=(a,b),1
            else: key,s=(b,a),-1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
T = pickle.load(open("field_tables.pkl","rb"))
T2 = T["T2"]
data = [("d",Fr(3,8),Fr(-1,4)), ("k",Fr(-3,8),Fr(1,4)),
        ("mr2h",Fr(-9,2),Fr(3)), ("r2h",Fr(-155,32),Fr(103,32)),
        ("x1",Fr(33,32),Fr(-21,32)), ("x2",Fr(-79,32),Fr(59,32))]
accP, accQ = {}, {}
for nm,p,q in data:
    vz = T2[nm][1]; vomz = T2["1-"+nm][1]
    w = wedge(vz,vomz)
    print(f"{nm}: vz={vz} vomz={vomz} wedge={w} c=({p},{q})")
    for (a,b),cab in w.items():
        for l3,c3 in vz.items():
            key=((a,b),l3)
            accP[key]=accP.get(key,Fr(0))+p*cab*c3
            accQ[key]=accQ.get(key,Fr(0))+q*cab*c3
            print(f"    (({a},{b}),{l3}): dP={p*cab*c3} dQ={q*cab*c3}")
print("P:", accP)
print("Q:", accQ)
# hand check for tensor r2:
handP = -Fr(3,8) + Fr(-9,2) + Fr(-155,32) + Fr(-79,32)
handQ = Fr(1,4) + 3 + Fr(103,32) + Fr(59,32)
print("hand ((r2,u),r2):", handP, handQ)
