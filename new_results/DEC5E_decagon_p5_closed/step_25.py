import sympy as sp
r5 = sp.sqrt(5)
s = (r5-1)/4                      # sin(pi/10)
c = sp.sqrt(10+2*r5)/4            # cos(pi/10)
m = 5
C5 = sp.expand(2880*m*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1))
C23 = sp.expand(-960*m*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1))
L = sp.sqrt((25-11*r5)/10)
def pq(e):
    e = sp.radsimp(sp.expand(e))
    p = sp.radsimp((e + e.subs(r5,-r5))/2); q = sp.radsimp((e - e.subs(r5,-r5))/(2*r5))
    return sp.simplify(p), sp.simplify(q)
# c/L in Q(sqrt5)?
print("c/L =", sp.nsimplify(sp.simplify(c/L), [r5]), "=", sp.simplify(c/L - (sp.Rational(5,2)+r5)))
v5 = sp.simplify(C5/L); v23 = sp.simplify(C23/L)
p5, q5 = pq(v5); p23, q23 = pq(v23)
print("C5/L  =", p5, "+", q5, "*sqrt5")
print("C23/L =", p23, "+", q23, "*sqrt5")
ok5 = (p5, q5) == (sp.Rational(73125,64), sp.Rational(-30825,64))
ok23 = (p23, q23) == (sp.Rational(-64875,64), sp.Rational(32775,64))
print("z5 coefficient matches DEC5D local block:", ok5)
print("z2z3 coefficient matches DEC5D local block:", ok23)
assert ok5 and ok23
open('progress.log','a').write("step_21: prop:order5 m=5 EXACT: C5^(5)/L = 73125/64-30825/64 sqrt5, C5^(2,3)/L = -64875/64+32775/64 sqrt5 = DEC5D local-block coefficients (sympy, exact)\n")
print("step_21 OK")
