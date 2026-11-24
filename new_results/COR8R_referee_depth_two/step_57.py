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
sum3 = S.S_add(S.symG((TWO,Z,B)), S.S_add(S.symG((Z,TWO,B)), S.symG((Z,B,TWO))))
Ep = expand(prod); Es = expand(sum3)
keys = sorted(set(Ep)|set(Es))
print("(e0=[2],e1=[c],e2=[1-c],e3=[1+c])")
print(f"{'component':15} {'prod':>6} {'sum3':>6} {'diff':>6}")
for k in keys:
    print(f"{str(k):15} {str(Ep.get(k,0)):>6} {str(Es.get(k,0)):>6} {str(Ep.get(k,0)+Es.get(k,0)):>6}")
print("\nNB: residual = prod + sum3 since expected prod = -sum3")
