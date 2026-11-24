import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc); Z = sp.Integer(0)
def expand(sym):
    out = {}
    for w,cf in sym.items():
        c0 = cf.get(0,Fraction(0))
        if not c0: continue
        partials = [((), c0)]
        for L in w:
            new = []
            opts = [(i, L[i]) for i in range(4) if L[i] != 0]
            for (pre, cpre) in partials:
                for (i, li) in opts:
                    new.append((pre+((i,),), cpre*li))
            partials = new
        for t,c in partials:
            key = tuple(x[0] for x in t)
            out[key] = out.get(key,Fraction(0))+c
    return {k:v for k,v in out.items() if v}
prod = S.S_shuffle(S.symG((TWO,)), S.symG((B,Z)))
parts = {'prod(G2*GB0)': prod, 's20B': S.symG((TWO,Z,B)), 's02B': S.symG((Z,TWO,B)), 's0B2': S.symG((Z,B,TWO))}
E = {k: expand(v) for k,v in parts.items()}
keys = sorted(set().union(*[set(e) for e in E.values()]))
print(f"{'comp':12}" + "".join(f"{k:>14}" for k in parts))
for kk in keys:
    print(f"{str(kk):12}" + "".join(f"{str(e.get(kk,0)):>14}" for e in E.values()))
print("\nexpect prod = -(s20B+s02B+s0B2): check prod + s20B + s02B + s0B2 = 0")
