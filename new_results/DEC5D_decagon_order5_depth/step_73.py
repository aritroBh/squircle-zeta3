
import sympy as sp
s5 = sp.sqrt(5)
phi = (1+s5)/2
s = 1/(2*phi)
cOverL = phi*(3*s5+5)/4     # cos(pi/10)/L  in Q(sqrt5)  [since L = sqrt(10-2 s5)*(3 s5-5)/20]
C5_5_over_L = sp.expand(2880*5*s**4*(65*s**6 - 63*s**4 + 17*s**2 - 1)*cOverL)
C5_23_over_L = sp.expand(-960*5*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)*cOverL)
def pq5(e):
    e = sp.expand(e)
    A = e.coeff(s5, 0); B = e.coeff(s5, 1)
    assert sp.expand(e - (A+B*s5)) == 0, e
    return (sp.Rational(A), sp.Rational(B))
p1 = pq5(C5_5_over_L); p2 = pq5(C5_23_over_L)
print("local-block zeta(5) coeff in P5/L:", p1)
print("local-block zeta(2)zeta(3) coeff in P5/L:", p2)
import mpmath as mp
mp.mp.dps = 40
print("numeric check:", mp.nstr(float(p1[0])+float(p1[1])*mp.sqrt(5), 25), "vs 65.5969467745153532529253385")
print("numeric check:", mp.nstr(float(p2[0])+float(p2[1])*mp.sqrt(5), 25), "vs 131.439499414931364059541672")
