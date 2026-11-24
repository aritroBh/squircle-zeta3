import cor8r_sym as S
import sympy as sp
cc = S.cc
half = sp.Rational(1,2)
# convB Li_{3,1}(1/2,1+c) = G(0,0,2,2/(1+c);1)
s31p_B = S.Li31_sym(half, 1+cc, 'B')
print("convB S(Li_{3,1}(1/2,1+c)) #terms:", len(s31p_B))
NAMES = ['-1','2','c','1-c','1+c']
def showL(L):
    parts=[]
    for n,e in zip(NAMES,L):
        if e: parts.append(f"{n}^{e}" if e!=1 else n)
    return "["+" ".join(parts)+"]"
for w,cf in sorted(s31p_B.items()):
    print("  ", " ".join(showL(L) for L in w), dict(cf))
