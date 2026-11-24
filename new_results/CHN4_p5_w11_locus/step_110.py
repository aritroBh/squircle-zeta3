import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
# Spot check P7.txt: reconstruct and verify irreducible, reversal-symmetric, and f11 relation at a test point
P7=sp.sympify(open('P7.txt').read())
print("P7 terms:", len(sp.Add.make_args(P7)), "total deg:", sp.Poly(P7,*A).total_degree())
print("P7 irreducible:", len(sp.factor_list(P7)[1])==1)
print("P7 reversal-symmetric:", sp.expand(P7 - P7.subs({a1:a5,a2:a4,a4:a2,a5:a1}))==0)
# verify f11 = -(1/10) a2a3a4 Sum P7 at a random point against engine
import sys; sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
def f11eng(a):
    a=[F(x) for x in a]; acc=F(0); T=PE.table(11)
    for pt in (a,a[::-1]):
        for w,pp in PE.volreg_words_point(pt,11)[11].items():
            c=T.reduce_word(w)[8]
            if c: acc+=pp*c
    return acc
import random
random.seed(42); ok=True
for _ in range(3):
    pt=tuple(random.randint(1,7) for _ in range(5))
    sa=sum(pt)
    expect=-F(1,10)*pt[1]*pt[2]*pt[3]*sa*F(P7.subs(dict(zip(A,pt))))
    got=f11eng(pt)
    match=(expect==got); ok=ok and match
    print(f"  pt {pt}: f11 engine={got}, -(1/10)a2a3a4(Sum)P7={expect}, match={match}")
print("P7.txt verified against engine:", ok)
