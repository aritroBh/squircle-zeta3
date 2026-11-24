from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, si, quad)
mp.dps = 100
def N(x,d=36): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
Lp = log(pi)
def bp3(n):  # d/ds b_n(s) at s=3
    return (2/pi)*( -pi**3*Lp*(-1)**n/n + 3*pi*(2*Lp+1)*(-1)**n/n**3 + 6*si(n*pi)/n**4 + 2*pi*(-1)**n/n**3 )
def bp3_quad(n): return (2/pi)*quad(lambda x: x**3*log(x)*sin(n*x), [0,pi])
print("bp3 formula check:", [nstr(abs(bp3(n)-bp3_quad(n)),3) for n in [1,2,4,6]])
gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
S4 = nsum(lambda k: gfun(2*k+1)*(-1)**k/(2*k+1)**4, [0, inf])
b4 = beta(4)
lhs = 6*b4 + (12/pi)*S4
rhs = pi**3*(2-log(2))/8
print("6 beta(4) + (12/pi) S4 =", N(lhs,40))
print("pi^3 (2 - log2)/8      =", N(rhs,40), " diff:", nstr(abs(lhs-rhs),3))
# hence S4 closed form:
S4closed = (pi/12)*(rhs - 6*b4)
print("S4 direct =", N(S4,40)); print("S4 closed =", N(S4closed,40), " diff:", nstr(abs(S4-S4closed),3))
# kink split at s=3: sum [bp3 - 6/n^4] sin(n pi/2) = T'(3) - 6 beta(4)
ks = nsum(lambda k: (bp3(2*k+1) - 6/(2*k+1)**4)*(-1)**k, [0, inf])
print("kink split T'(3):", N(ks,40), " T'(3)-6beta(4) =", N((pi/2)**3*log(pi/2)-6*b4,40), " diff:", nstr(abs(ks-((pi/2)**3*log(pi/2)-6*b4)),3))

print()
print("== kink amplitude census: n^{s+1} b_n(s) -> (2/pi) Gamma(s+1) cos(pi s/2) ;  n^{s+1} a_n(s) -> -(2/pi) Gamma(s+1) sin(pi s/2) ==")
def a_g(n,s): return (2/pi)*quad(lambda x: x**s*cos(n*x), [0,pi])
def b_g(n,s): return (2/pi)*quad(lambda x: x**s*sin(n*x), [0,pi])
nn = 400
for s in [mpf(2), mpf('2.5'), mpf(3), mpf('3.5'), mpf(4)]:
    ka = -(2/pi)*gamma(s+1)*sin(pi*s/2)
    kb = (2/pi)*gamma(s+1)*cos(pi*s/2)
    va = a_g(nn,s)*nn**(s+1); vb = b_g(nn,s)*nn**(s+1)
    print(f"s={nstr(s,3)}: a-amp pred {nstr(ka,10)}, n^{{s+1}} a_n({nn}) = {nstr(va,10)} | b-amp pred {nstr(kb,10)}, n^{{s+1}} b_n = {nstr(vb,10)}")
logline("step_38: T'(3)/beta(4)/S4 identity verified; kink amplitude census table")
