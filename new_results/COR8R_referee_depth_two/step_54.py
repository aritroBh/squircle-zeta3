import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc); Z = sp.Integer(0)
s2   = S.symG((TWO,))
sB0  = S.symG((B,Z))
s20B = S.symG((TWO,Z,B))
s02B = S.symG((Z,TWO,B))
s0B2 = S.symG((Z,B,TWO))
prod = S.S_shuffle(s2, sB0)
res = S.S_add(prod, S.S_add(s20B, S.S_add(s02B, s0B2)))
# expand multilinearly: each letter is a vector in Z^4; tensor = multilinear form
exp = {}
for w,cf in res.items():
    c0 = cf.get(0,Fraction(0))
    if not c0: continue
    terms = [((), c0)]
    for L in w:
        new = []
        for idx,v in enumerate(L):
            if v:
                for t,c in terms:
                    new.append((t+( (idx,v) ),c))
        # need full multilinear expansion: distribute each letter
        pass
    # proper distribution:
    from itertools import product as iproduct
    choices = []
    for L in w:
        opts = [(i,)*1 for i in range(4) for _ in [0]]
        opts = []
        for i in range(4):
            if L[i]: opts.append((i, L[i]))
        choices.append(opts if opts else [None])
    def rec(k, acc, coeff):
        if k==len(w):
            exp[tuple(acc)] = exp.get(tuple(acc),Fraction(0))+coeff
            return
        if choices[k] is [None] or choices[k]==[None]:
            rec(k+1, acc, coeff); return
        for (i,vi) in choices[k]:
            # letter = sum_i L[i] e_i, so choose coordinate i with multiplicity L[i]
            for _ in range(abs(vi)):
                pass
            acc2 = acc+[(i,)]*abs(vi)
            rec(k+1, acc2, coeff*(1 if vi>0 else (-1)**abs(vi)))
    rec(0, [], c0)
resid = {k:v for k,v in exp.items() if v}
print("multilinear expansion nonzero components:", len(resid))
for k,v in sorted(resid.items()):
    print("  ", k, v)
