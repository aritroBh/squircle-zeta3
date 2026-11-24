code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, quad, fsum, inf, pi
setdps(80)
s = 4
def Jdirect(u):
    tot = mpf(0); d = 0
    while True:
        d += 1
        t = log(1+exp(-d/u))
        tot += t
        if t < mpf(10)**(-80): break
    return tot
def Jmod(u):
    # pi^2 u/12 - ln2/2 + 1/(24u) - sum_n ln(1+e^{-2 pi^2 u n})
    tot = pi**2*u/12 - log(2)/2 + 1/(24*u)
    n = 0
    while True:
        n += 1
        t = log(1+exp(-2*pi**2*u*n))
        tot -= t
        if t < mpf(10)**(-80): break
    return tot
# consistency check
for u in [0.7, 1.5, 3.0]:
    print(u, nstr(Jdirect(u), 20), nstr(Jmod(u), 20))
# continuum integral: [0,2] direct; [2,inf): analytic asymptote + quad remainder
I1 = quad(lambda u: Jdirect(u)/u**s, [0, 2])
def rem(u):
    return -(fsum([log(1+exp(-2*pi**2*u*n)) for n in range(1,20)]))/u**s
I2a = quad(rem, [2, inf])
# analytic: int_2^inf u^{-s}(pi^2 u/12 - ln2/2 + 1/(24u)) du
I2b = pi**2/12 * power(2, 2-s)/(s-2) - log(2)/2 * power(2, 1-s)/(s-1) + mpf(1)/24 * power(2, -s)/s
print("I1:", nstr(I1, 25))
print("I2a:", nstr(I2a, 25), " I2b:", nstr(I2b, 25))
Itot = I1 + I2a + I2b
print("continuum int:", nstr(Itot, 30))
pred = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
print("modular pred   :", nstr(pred, 30))
print("fit value      : 2.2881505889")
'''
exec(code)
