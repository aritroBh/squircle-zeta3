# Check: P(p) = 4*sqrt2*(1+eta) * int_0^1 (t(1-t))^(eta/2) sqrt(cosh(eta*log(t/(1-t)))) dt, eta = 1/p - 1,
# and the Taylor coefficients predicted by the Beta/Gamma operator formula
#   P = sqrt2 * A(mu) * sum_k c_k (2mu)^(2k) [d_s^(2k) Gamma(1+mu+s)Gamma(1+mu-s)]_{s=0} / Gamma(1+mu)^2,   mu = eta/2
from mpmath import mp, mpf, quad, sqrt, cosh, log, gamma, pi, zeta, taylor, diff
mp.dps = 50
def P_direct(p):
    e = 1/mpf(p); f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
    return 8*quad(f, [0, mpf(1)/8, mpf(1)/4, mpf(3)/8, mpf(1)/2])
# Taylor coefficients of sqrt(cosh z) = sum c_k z^(2k)
cz = taylor(lambda z: sqrt(cosh(z)), 0, 12)
def P_formula(p, K=6):
    eta = 1/mpf(p) - 1; mu = eta/2
    G = lambda s: gamma(1+mu+s)*gamma(1+mu-s)
    tot = 0
    for k in range(K+1):
        dk = G(0) if k == 0 else diff(G, 0, 2*k)
        tot += cz[2*k]*eta**(2*k)*dk
    return 4*sqrt(2)*(1+eta)*tot/gamma(2+eta)
for p in ['1.01', '0.99', '1.05', '1.1']:
    pd, pf = P_direct(p), P_formula(p)
    print(p, mp.nstr(pd, 25), " formula-direct =", mp.nstr(pf-pd, 3))
# Closed-form Taylor coefficients in d = p-1 (from the formula) vs. claim
d = mpf(10)**-3
print("P(1+d) - [4sqrt2 + (sqrt2/6)pi^2 d^2 + (sqrt2 zeta3 - sqrt2 pi^2/3) d^3] =", mp.nstr(P_direct(1+d) - (4*sqrt(2) + sqrt(2)/6*pi**2*d**2 + (sqrt(2)*zeta(3) - sqrt(2)*pi**2/3)*d**3), 5), " (expect ~d^4 = 1e-12)")
