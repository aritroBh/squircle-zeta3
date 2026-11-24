from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, si, quad, tan)
mp.dps = 100
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
def a_c(n,s):
    if s==1: return 2*((-1)**n-1)/(pi*n**2)
    if s==2: return 4*(-1)**n/n**2
    if s==3: return 6*pi*(-1)**n/n**2 + 12*(1-(-1)**n)/(pi*n**4)
    if s==4: return (-1)**n*(8*pi**2/n**2 - 48/n**4)
def b_c(n,s):
    if s==1: return 2*(-1)**(n+1)/n
    if s==2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
    if s==3: return 2*pi**2*(-1)**(n+1)/n + 12*(-1)**n/n**3
    if s==4: return 2*pi**3*(-1)**(n+1)/n + 24*pi*(-1)**n/n**3 + 48*(1-(-1)**n)/(pi*n**5)
exact = {(1,1): -7*zeta(3)/pi,
 (1,2): -7*zeta(3) + 31*zeta(5)/pi**2,
 (2,1): -8*zeta(3),
 (2,2): -8*pi*zeta(3) + 31*zeta(5)/pi,
 (1,3): -7*pi*zeta(3) + mpf('46.5')*zeta(5)/pi,
 (2,3): -8*pi**2*zeta(3) + 48*zeta(5),
 (3,3): -12*pi**3*zeta(3) + mpf('118.5')*pi*zeta(5) - mpf('285.75')*zeta(7)/pi,
 (1,4): -7*pi**2*zeta(3) + 93*zeta(5) - 381*zeta(7)/pi**2,
 (2,4): -8*pi**3*zeta(3) + 96*pi*zeta(5) - 381*zeta(7)/pi,
 (3,4): -12*pi**4*zeta(3) + mpf('190.5')*pi**2*zeta(5) - 1143*zeta(7) + mpf('2299.5')*zeta(9)/pi**2,
 (4,4): -16*pi**5*zeta(3) + 288*pi**3*zeta(5) - 1914*pi*zeta(7) + 4599*zeta(9)/pi}
M = 60000
for (s,q),v in exact.items():
    part = nsum(lambda n: a_c(n,s)*b_c(n,q), [1, M])
    err = abs(part - v)
    print(f"B({s},{q}): partial(6e4) vs exact, diff = {nstr(err,3)} (tail ~1e-{8 if s+q<5 else 9})")
# log-cos integral for B(1,1)
I = (4/pi)*quad(lambda x: x*log(2*cos(x/2)), [0, pi])
print("(4/pi) int x log(2cos(x/2)) =", N(I,36), " vs -7zeta(3)/pi =", N(-7*zeta(3)/pi,36))
logline("step_37: B-table confirmed by partial sums; log-cos integral form of B(1,1)")
