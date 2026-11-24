from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, si, euler, pslq, quad, tan)
mp.dps = 110
def N(x,d=36): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi)
# R4 = sum (Si(n pi) - pi/2)/n^4 : direct nsum (smooth alternating-free? (Si-pi/2) oscillates: use g-reindex pieces)
# (Si(nπ)-π/2) = (−1)^{n+1} g(n), g smooth: so R4 = sum (-1)^{n+1} g(n)/n^4 : alternating smooth
gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
R4 = nsum(lambda n: (-1)**(n+1)*gfun(n)/n**4, [1, inf])
print("R4 =", N(R4,50))
def bp3(n):  return (2/pi)*( -pi**3*Lp*(-1)**n/n + 3*pi*(2*Lp+1)*(-1)**n/n**3 + 6*si(n*pi)/n**4 + 2*pi*(-1)**n/n**3 )
Wp3_pieces = 2*pi**2*Lp*log(2) - (9*Lp+mpf('7.5'))*zeta(3) + pi**4/15 + (12/pi)*R4
Wp3_direct = nsum(lambda n: bp3(n), [1, inf])
Wp3_int = (1/pi)*quad(lambda t: t**3*log(t)/tan(t/2), [0, pi])
print("W'(3) pieces =", N(Wp3_pieces,45))
print("W'(3) direct =", N(Wp3_direct,45))
print("W'(3) integral =", N(Wp3_int,45))
# PSLQ R4 against plausible basis
basis = [mpf(1), log(2), Lp, euler, pi**2, pi**4, zeta(3), zeta(3)*Lp, euler*zeta(3), zeta(3,derivative=1), zeta(5), pi**4*log(2), pi**4*Lp, zeta(4,derivative=1)]
names = ['1','log2','Lp','gam','pi^2','pi^4','z3','z3 Lp','gam z3',"z3'",'z5','pi^4 log2','pi^4 Lp',"z4'"]
rel = pslq([R4]+basis, tol=mpf('1e-95'), maxcoeff=10**10, maxsteps=10000)
print("PSLQ R4:", rel)
if rel:
    resid = abs(sum(c*v for c,v in zip(rel,[R4]+basis)))
    print("  residual:", nstr(resid,3))
logline("step_49: R4 and W'(3) computed; PSLQ R4 null expected -> OPEN")
