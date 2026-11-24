
import sympy as sp
s5 = sp.sqrt(5)
phi = (1+s5)/2
s = 1/(2*phi)                  # sin(pi/10)
# c/L factor: cos(pi/10) = phi*sqrt(10-2 s5)/4 ; L = sqrt((25-11 s5)/10)
ccos = phi*sp.sqrt(10-2*s5)/4
L = sp.sqrt((25-11*s5)/10)
C5_5 = sp.expand(2880*5*s**4*ccos*(65*s**6 - 63*s**4 + 17*s**2 - 1))
C5_23 = sp.expand(-960*5*ccos*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1))
r5_ = sp.Symbol
# C/L: the radical sqrt(10-2 s5) cancels against L? L is a different radical; expand and simplify
import sympy as sy
q1 = sp.radsimp(sp.expand(C5_5/L))
q2 = sp.radsimp(sp.expand(C5_23/L))
# extract (p,q) with p+q sqrt5
def pq5(e):
    e = sp.expand(e)
    A = e.coeff(s5, 0); B = e.coeff(s5, 1)
    assert sp.expand(e - (A+B*s5)) == 0, e
    return (sp.Rational(A), sp.Rational(B))
print("zeta(5) coeff of local block in P5/L  =", pq5(q1))
print("zeta(2)zeta(3) coeff of local block in P5/L =", pq5(q2))
import mpmath as mp
mp.mp.dps = 40
print("numeric:", mp.nstr(float(pq5(q1)[0])+float(pq5(q1)[1])*mp.sqrt(5), 30), "vs", mp.nstr(65.59694677451535325292533845872515691473, 30))
print("numeric:", mp.nstr(float(pq5(q2)[0])+float(pq5(q2)[1])*mp.sqrt(5), 30), "vs", mp.nstr(131.4394994149313640595416717604309158838, 30))
