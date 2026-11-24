import cor8r_sym as S
import sympy as sp
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc); Z = sp.Integer(0)
s2   = S.symG((TWO,))
sB0  = S.symG((B,Z))
s20B = S.symG((TWO,Z,B))
s02B = S.symG((Z,TWO,B))
s0B2 = S.symG((Z,B,TWO))
prod = S.S_shuffle(s2, sB0)
res = S.S_add(prod, S.S_add(s20B, S.S_add(s02B, s0B2)))
NAMES = ['2','c','1-c','1+c']
def showL(L):
    parts=[]
    for n,e in zip(NAMES,L):
        if e: parts.append(f"{n}^{e}" if e!=1 else n)
    return "["+" ".join(parts)+"]"
print("residual #terms:", len(res))
for w,cf in res.items():
    print("  ", " ".join(showL(L) for L in w), dict(cf))
print()
print("s2 =", {w:dict(cf) for w,cf in s2.items()})
print("sB0 =", {w:dict(cf) for w,cf in sB0.items()})
print("prod =", {w:dict(cf) for w,cf in prod.items()})
