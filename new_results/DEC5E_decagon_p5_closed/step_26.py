import sympy as sp
r5 = sp.sqrt(5)
s = (r5-1)/4
m = 5
# C5 = A5 * c, C23 = A23 * c with A5, A23 in Q(sqrt5)
A5 = sp.expand(2880*m*s**4*(65*s**6 - 63*s**4 + 17*s**2 - 1))
A23 = sp.expand(-960*m*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1))
print("A5 =", A5, "   A23 =", A23)
# c/L = 5/2+sqrt5: prove by squares
c = sp.sqrt(10+2*r5)/4; L = sp.sqrt((25-11*r5)/10)
lhs2 = sp.radsimp(sp.expand(c**2/L**2)); rhs2 = sp.expand((sp.Rational(5,2)+r5)**2)
print("(c/L)^2 =", lhs2, "= (5/2+sqrt5)^2 =", rhs2, "equal:", sp.simplify(lhs2-rhs2)==0, "; both positive")
def pq(e):
    e = sp.expand(e)
    p = sp.expand((e + e.subs(r5,-r5))/2); q = sp.expand((e - e.subs(r5,-r5))/(2*r5))
    return sp.simplify(p), sp.simplify(q)
v5 = sp.expand(A5*(sp.Rational(5,2)+r5)); v23 = sp.expand(A23*(sp.Rational(5,2)+r5))
p5,q5 = pq(v5); p23,q23 = pq(v23)
print("C5^(5)/L  =", p5, "+", q5, "sqrt5")
print("C5^(2,3)/L =", p23, "+", q23, "sqrt5")
assert (p5,q5)==(sp.Rational(73125,64), sp.Rational(-30825,64))
assert (p23,q23)==(sp.Rational(-64875,64), sp.Rational(32775,64))
print("EXACT MATCH with DEC5D local-block z5, z2z3 coefficients")
open('progress.log','a').write("step_26: prop:order5 m=5 EXACT (radical-free): C5^(5)/L=73125/64-30825/64*sqrt5, C5^(2,3)/L=-64875/64+32775/64*sqrt5 = DEC5D local-block coefficients\n")
