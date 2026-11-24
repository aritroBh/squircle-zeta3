import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, gamma, sin, pi
def F(y):
    return quad(lambda w: w**(-y) * ((1+w)/sqrt(1+w**2) - 1)/w, [0,1])
def E(y):
    return gamma(1+y/2)*gamma(mpf('0.5')-y/2)/sqrt(pi)
for y in [mpf('0.3'), mpf('0.13'), mpf('0.71')]:
    lhs = (F(y)+F(-y))/2
    rhs = (E(y)-E(-y))/(2*y)
    print(y, mpmath.nstr(lhs,25), mpmath.nstr(rhs,25), mpmath.nstr(lhs-rhs,8))
# also check the odd part: (F(y)-F(-y))/2 vs (E(y)+E(-y))/(2y) - 1/y?
for y in [mpf('0.3'), mpf('0.71')]:
    lhs = (F(y)-F(-y))/2
    rhs = (E(y)+E(-y))/(2*y) - 1/y   # even part of E minus e_0
    rhs2 = (E(y)+E(-y))/(2*y)
    print('odd', y, mpmath.nstr(lhs,25), mpmath.nstr(rhs,25), mpmath.nstr(rhs2,25))
