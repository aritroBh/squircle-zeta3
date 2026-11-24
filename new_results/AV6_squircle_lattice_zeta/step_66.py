code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, quad, fsum, inf
setdps(80)
s = 4
def J(u):
    # sum_{d>=1} ln(1+e^{-d/u})
    if u <= 0: return mpf(0)
    tot = mpf(0); d = 0
    while True:
        d += 1
        t = log(1+exp(-d/u))
        tot += t
        if t < mpf(10)**(-90): break
    return tot
I = quad(lambda u: J(u)/u**s, [0, 1, inf])
print("continuum integral:", nstr(I, 30))
pred = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
print("modular prediction :", nstr(pred, 30))
print("fit value          : 2.2881505889...")
# check the Mellin evaluation: int v^{s-2} ln(1+e^{-v}) dv = Gamma(s-1) eta(s)
M = quad(lambda v: v**(s-2)*log(1+exp(-v)), [0, inf])
print("L_1(s-1) =", nstr(M, 30), " vs Gamma eta:", nstr(gamma(s-1)*(1-power(2,-(s-1)))*zeta(s), 30))
'''
exec(code)
