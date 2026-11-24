import sympy as sp
from fractions import Fraction as F

a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
A = [a1,a2,a3,a4,a5]
SA = a1+a2+a3+a4+a5

def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator * sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))

P8s = dict_to_sym(P8)
C4 = sp.expand(10*P8s / (a2*a3*a4*SA))
assert sp.expand(a2*a3*a4*SA*C4/10 - P8s) == 0
print("C4: degree", sp.Poly(C4, a1,a2,a3,a4,a5).total_degree(), "terms:", len(sp.Poly(C4).terms()))
print("C4(1,1,1,1,1) =", C4.subs({a1:1,a2:1,a3:1,a4:1,a5:1}))

P10s = dict_to_sym(P10_37)
Q7p = sp.expand(224*P10s / (a2*a3*a4*SA))
assert sp.expand(a2*a3*a4*SA*Q7p/224 - P10s) == 0
print("Q7': degree", sp.Poly(Q7p).total_degree(), "terms:", len(sp.Poly(Q7p).terms()))

# also check z2z(3,5) form: -a2a3a4(SA)(a2+a4)(a1+a3+a5)C4/10
P235s = dict_to_sym(P10_235)
check = sp.expand(P235s + a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10)
print("z2z(3,5) identity check:", check == 0)
