import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc); Z = sp.Integer(0)
syms = {
 'prod': S.S_shuffle(S.symG((TWO,)), S.symG((B,Z))),
 's20B': S.symG((TWO,Z,B)),
 's02B': S.symG((Z,TWO,B)),
 's0B2': S.symG((Z,B,TWO)),
}
res = S.S_add(syms['prod'], S.S_add(syms['s20B'], S.S_add(syms['s02B'], syms['s0B2'])))
def expand(sym):
    out = {}
    for w,cf in sym.items():
        c0 = cf.get(0,Fraction(0))
        if not c0: continue
        partials = [((), c0)]
        for L in w:
            new = []
            coords = []
            for i in range(4):
                coords += [i]*abs(L[i])
            sign = 1 if all(L[i]>=0 for i in range(4)) else None
            # letter = sum L[i] e_i; multilinear: tensor factor expands to sum over choices with multiplicities L[i]
            # L[i] can be negative: e_i^(-1) means coefficient -1 with one e_i? NO:
            # [f^(-1)] = -[f]; so letter vector L corresponds to sum_i L[i]*e_i as a VECTOR.
            # multilinear expansion of tensor with vector entries: entry v expands as sum_i L[i]*(e_i)
            opts = []
            for i in range(4):
                if L[i] != 0:
                    opts.append((i, L[i]))
            for (pre, cpre) in partials:
                for (i, li) in opts:
                    new.append((pre+((i,),), cpre*li))
            partials = new
        for t,c in partials:
            key = tuple(x[0] for x in t)
            out[key] = out.get(key,Fraction(0))+c
    return {k:v for k,v in out.items() if v}
E = expand(res)
print("residual expansion:", len(E))
for k,v in sorted(E.items()): print("  ", k, v)
# Also expand each piece separately to see the mismatch structure
for name,s in syms.items():
    E2 = expand(s)
    print(f"\n{name}: {len(E2)} components")
