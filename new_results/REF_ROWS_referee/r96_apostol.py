# R96: Apostol's integral over the quarter squircle.
from mpmath import mp, mpf, gamma, loggamma, zeta, pi, log, sqrt, quad, nsum, inf, nstr, exp, polylog, catalan
mp.dps = 40
def Zsum(p):
    return nsum(lambda m: exp(2*loggamma(1+m/p)-loggamma(1+2*m/p))/m**2, [1, inf])
def Zquad(p):
    # inner integral over y in closed form: int_0^Y dy/(1-xy) = -log(1-xY)/x
    p = mpf(p)
    def f(x):
        Y = (1-x**p)**(1/p)
        return -log(1-x*Y)/x
    return quad(f, [0, mpf(1)/2**(1/p), 1])
phi = (1+sqrt(5))/2
for p, closed, name in [(mpf(1)/2, pi**2/9-4*log(phi)**2, "pi^2/9-4log^2 phi"), (mpf(1), zeta(2)/3, "zeta(2)/3"), (mpf(2), 7*zeta(2)/12, "7zeta(2)/12"), (mpf(3), None, ""), (mpf('7.5'), None, "")]:
    s = Zsum(p); q = Zquad(p)
    line = "p=%s  sum=%s  |sum-quad|=%s" % (nstr(p, 4), nstr(s, 30), nstr(abs(s-q), 3))
    if closed is not None: line += "  |sum-%s|=%s" % (name, nstr(abs(s-closed), 3))
    print(line)
# epsilon expansion: G(x)=Gamma(1+x)^2/Gamma(1+2x)
# epsilon expansion: see r96_integral.py (series for 1-G near 0 avoids cancellation)
