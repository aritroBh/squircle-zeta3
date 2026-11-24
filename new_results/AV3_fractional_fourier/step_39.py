from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta)
mp.dps = 80
def N(x,d=42): return nstr(x,d)
def a_c(n,s):
    if s==1: return 2*((-1)**n-1)/(pi*n**2)
    if s==2: return 4*(-1)**n/n**2
def b_c(n,s):
    if s==1: return 2*(-1)**(n+1)/n
    if s==2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
M = 200000
p11 = nsum(lambda n: a_c(n,1)*b_c(n,1), [1, M])
p12 = nsum(lambda n: a_c(n,1)*b_c(n,2), [1, M])
print("partial B(1,1) to 2e5:", N(p11,30), " tail ~ O(1/M^2) ~ 1e-11")
print("partial B(1,2) to 2e5:", N(p12,30))
print("-7 zeta(3)/pi          :", N(-7*zeta(3)/pi,30))
print("-(4/pi)(zeta(3)+eta(3)):", N(-(4/pi)*(zeta(3)+mpf(3)/4*zeta(3)),30))
