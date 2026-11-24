# Corrected: 1-G(x) computed without cancellation, via the series log G = sum_{k>=2} (-1)^k zeta(k)(2-2^k)/k x^k
# for small x, and expm1 elsewhere.
from mpmath import mp, mpf, loggamma, exp, expm1, quad, inf, nstr, zeta, nsum, fsum
mp.dps = 50
def logG(x):
    if x < mpf('0.05'):
        return fsum((-1)**k*zeta(k)*(2-2**k)/k*x**k for k in range(2, 90))
    return 2*loggamma(1+x)-loggamma(1+2*x)
f = lambda x: -expm1(logG(x))/x**2
v, err = quad(f, [0, mpf('0.05'), mpf(1)/2, 1, 2, 4, 8, 16, 32, 64, 128, inf], error=True)
print("int_0^inf (1-G)/x^2 =", nstr(v, 40), " err est", nstr(err, 3), "  printed 1.7330025902")
mp.dps = 40
for p in [50, 100, 200, 400]:
    e = mpf(1)/p
    Z = nsum(lambda m: exp(2*loggamma(1+m*e)-loggamma(1+2*m*e))/m**2, [1, inf])
    r = (zeta(2)-Z - e*v + zeta(2)/2*e**2)/e**3
    print("p=%d  (zeta2-Z-eI+zeta2 e^2/2)/e^3 = %s ; zeta3/6=%s ; (r-zeta3/6)*p^2 = %s" % (p, nstr(r, 10), nstr(zeta(3)/6, 10), nstr((r-zeta(3)/6)*p**2, 8)))
from mpmath import taylor
g5 = taylor(lambda x: exp(logG(x)), 0, 5)[5]
print("predicted limit of (r-zeta3/6)*p^2 = -g5*zeta(-3) =", nstr(-g5*zeta(-3), 10))
