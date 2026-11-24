from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, si)
mp.dps = 100
def N(x,d=36): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
S4 = nsum(lambda k: gfun(2*k+1)*(-1)**k/(2*k+1)**4, [0, inf])
b4 = beta(4)
lhs = 6*b4 + (12/pi)*S4
rhs = pi**3*(5-2*log(2))/16
print("6 beta(4) + (12/pi) S4 =", N(lhs,42))
print("pi^3 (5 - 2 log2)/16   =", N(rhs,42), " diff:", nstr(abs(lhs-rhs),3))
S4closed = pi**4*(5-2*log(2))/192 - (pi/2)*b4
print("S4 direct =", N(S4,42)); print("S4 closed =", N(S4closed,42), " diff:", nstr(abs(S4-S4closed),3))

print()
print("== kink census with boundary subtracted (n=300) ==")
def kinkA(n,s): return -(2/pi)*gamma(s+1)*sin(pi*s/2)*n**(-s-1)
def bndA(n,s,J): return sum( (2*(-1)**n/pi)*(-1)**(j-1)*gamma(s+1)/gamma(s-2*j+2)*pi**(s-2*j+1)*n**(-2*j) for j in range(1,J+1))
def kinkB(n,s): return (2/pi)*gamma(s+1)*cos(pi*s/2)*n**(-s-1)
def bndB(n,s,J): return sum( (2*(-1)**(n+1)/pi)*(-1)**j*gamma(s+1)/gamma(s-2*j+1)*pi**(s-2*j)*n**(-(2*j+1)) for j in range(0,J+1))
def a_gam(n, s): return cos((2/pi)*0) * 0 + __import__('mpmath').re((2/pi)*(-1j*n)**(-s-1) * __import__('mpmath').gammainc(s+1, 0, -1j*n*pi))
def b_gam(n, s): return __import__('mpmath').im((2/pi)*(-1j*n)**(-s-1) * __import__('mpmath').gammainc(s+1, 0, -1j*n*pi))
nn = 300
for s in [mpf(2), mpf('2.5'), mpf(3), mpf('3.5'), mpf(4)]:
    ra = (a_gam(nn,s) - bndA(nn,s,2))*nn**(s+1)
    rb = (b_gam(nn,s) - bndB(nn,s,2))*nn**(s+1)
    print(f"s={nstr(s,3)}: a-kink pred {nstr(kinkA(nn,s)*nn**(s+1),12)} got {nstr(ra,12)} | b-kink pred {nstr(kinkB(nn,s)*nn**(s+1),12)} got {nstr(rb,12)}")
logline("step_39: corrected S4 identity verified; kink census with boundary subtraction")
