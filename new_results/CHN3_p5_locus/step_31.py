exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
f = sp.factor_list(R, a1,a2,a4)
G = None; lin = []
for g,m in f[1]:
    if sp.Poly(g,a1,a2,a4).total_degree() > 2: G = sp.expand(g)
    else: lin.append(g)
print("linear factors:", lin, " scalar:", f[0])
print("G degree:", sp.Poly(G,a1,a2,a4).total_degree())
# Rehomogenize: check that a1*(a1+a3)*G_hom has degree 24 and matches homogeneous resultant at random point
# G_hom: homogenize with a3
Gh = sp.expand(sum(c * a1**e[0]*a2**e[1]*a4**e[2] * a3**(22-sum(e)) for e,c in sp.Poly(G,a1,a2,a4).terms()))
# verify: evaluate homogeneous resultant directly at a random point via roots
import random
random.seed(1)
pt = {a1:sp.Rational(3,2), a2:sp.Rational(2,1), a3:sp.Rational(1,1), a4:sp.Rational(5,2)}
C4,Q7 = get_C4_Q7()
import sympy.polys.domains as D
Rh = sp.resultant(C4, Q7, a5)  # expensive?? NO - this is the full 4-var resultant, skip
print("skip full homogeneous resultant (expensive); instead verify factor values:")
val = sp.expand(R - 16*a1*(a1+1)*G)
print("R == 16*a1*(a1+1)*G dehomogenized:", val==0)
