import cor8r_sym as S
import sympy as sp
cc = S.cc
A = 1/(1+cc); B = 2/(1+cc)
NAMES = ['-1','2','c','1-c','1+c']
def showL(L):
    parts=[]
    for n,e in zip(NAMES,L):
        if e: parts.append(f"{n}^{e}" if e!=1 else n)
    return "["+" ".join(parts)+"]"
def show(s):
    for w,cf in sorted(s.items()):
        print("   ", " ".join(showL(L) for L in w), dict(cf))
Z = sp.Integer(0)
print("symG((A,0,B)):"); show(S.symG((A,Z,B)))
print("manual: -[1-c-2]... expect -[(c-1)/2] [B] [c]: letters (1,-1,0,1,0),(0,1,0,0,-1),(0,0,1,0,0), coeff -1")
print()
print("symG((0,A,B)):"); show(S.symG((Z,A,B)))
print("symG((0,B,A)):"); show(S.symG((Z,B,A)))
print("symG((B,0,A)):"); show(S.symG((B,Z,A)))
print("symG((A,B,A)):"); show(S.symG((A,B,A)))
print("symG((B,A,B)):"); show(S.symG((B,A,B)))
print("symG((A,B)):"); show(S.symG((A,B)))
print("symG((B,A)):"); show(S.symG((B,A)))
print("symG((0,A)):"); show(S.symG((Z,A)))
print("symG((0,B)):"); show(S.symG((Z,B)))
