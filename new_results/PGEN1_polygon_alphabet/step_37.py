# step_14: delta(N3) for the generic octagon, edge by edge, two ways (OCT9d basis vs CR basis)
import os, sys, pickle
from fractions import Fraction as Fr
import sympy as sp
base = "<project>"

def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def vmul(v,a):
    a=Fr(a); return {k:c*a for k,c in v.items() if c*a}
def vecQ(q):
    q=Fr(q); out={}
    for m,sgn in ((abs(q.numerator),1),(q.denominator,-1)):
        for p,ex in sp.factorint(int(m)).items():
            out[p]=out.get(p,Fr(0))+sgn*ex
    return {k:c for k,c in out.items() if c}
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            if a<b: key,s=(a,b),1
            else: key,s=(b,a),-1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(terms):
    out={}
    for c,z in terms:
        c=Fr(c); vz,vomz=vecQ(z),vecQ(1-z)
        w=wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                out[key]=out.get(key,Fr(0))+c*cab*c3
    return {k:v for k,v in out.items() if v}

sf = pickle.load(open(os.path.join(base,"research/08_theory/OCT9d_octagon_short/short_form.pkl"),"rb"))
G8 = [Fr(s) for s in sf['G8']]
B0 = [Fr(s) for s in sf['B0']]
oct9d = {}
for e in (1,2,3):
    coeffs = [Fr(s) for s in sf['sols123'][e][:8]]
    oct9d[e] = [(c,z) for c,z in zip(coeffs,G8) if c]
coeffs0 = [Fr(s) for s in sf['sol0'][:20]]
oct9d[0] = [(c,z) for c,z in zip(coeffs0,B0) if c]

CRp = pickle.load(open("CR_presentations.pkl","rb"))

deltas = {}
for e in range(4):
    d1 = delta_li3(oct9d[e])
    d2 = delta_li3(CRp[e]['presentation'])
    same = vadd(d1, vmul(d2,-1)) == {}
    deltas[e] = d1
    print(f"=== edge {e}: delta from OCT9d basis ({len(oct9d[e])} Li3) vs CR basis ({len(CRp[e]['presentation'])} Li3): identical = {same}")
    if not same:
        print("  DIFF:", vadd(d1, vmul(d2,-1)))
    for k,v in sorted(d1.items()):
        print(f"   {v} * ({k[0][0]}^{k[0][1]}) (x) {k[1]}")
tot = vadd(*deltas.values())
print("=== TOTAL delta(N3) over the four edges (length factors are Q-independent scalars):")
print("  nonzero per edge:", {e: bool(deltas[e]) for e in range(4)})
print("  sum (formal, ignoring independent scalar factors):", tot if tot else "0")
pickle.dump(deltas, open("delta_oct9d.pkl","wb"))
