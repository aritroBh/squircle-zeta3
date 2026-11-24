import cor8r_sym as S
import sympy as sp
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc)
Z = sp.Integer(0)
NAMES = ['-1','2','c','1-c','1+c']
def showL(L):
    parts=[]
    for n,e in zip(NAMES,L):
        if e: parts.append(f"{n}^{e}" if e!=1 else n)
    return "["+" ".join(parts)+"]"
s2  = S.symG((TWO,))
sB  = S.symG((B,))
sB0 = S.symG((B,Z))     # L=2 solved
s20B = S.symG((TWO,Z,B))
s02B = S.symG((Z,TWO,B))
s0B2 = S.symG((Z,B,TWO))
lhs = S.S_shuffle(s2, sB0)
rhs = S.S_add(S.S_add(s20B, s02B), s0B2)
res = S.S_add(lhs, rhs)   # expect 0 if G(2)G(B,0) = -(G(2,0,B)+G(0,2,B)+G(0,B,2))
print("S(G(2)):"); [print("  ", " ".join(showL(L) for L in w), cf) for w,cf in s2.items()]
print("S(G(B,0)):"); [print("  ", " ".join(showL(L) for L in w), cf) for w,cf in sB0.items()]
print("residual #terms:", len(res))
for w,cf in list(res.items())[:20]:
    print("  ", " ".join(showL(L) for L in w), cf)
