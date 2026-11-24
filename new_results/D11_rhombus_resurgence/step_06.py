# Replicate rhombus_formula.py checks inline
from mpmath import mp, mpf, quad, sqrt, cosh, log, gamma, pi, zeta, taylor, diff
mp.dps = 50
def P_direct(p):
    e = 1/mpf(p); f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
    return 8*quad(f, [0, mpf(1)/8, mpf(1)/4, mpf(3)/8, mpf(1)/2])
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
d = mpf(10)**-3
print("resid d^3 check:", mp.nstr(P_direct(1+d) - (4*sqrt(2) + sqrt(2)/6*pi**2*d**2 + (sqrt(2)*zeta(3) - sqrt(2)*pi**2/3)*d**3), 5))

# Also check the l-variable weight w_eta claim: F(eta) = int_R e^{-l(1+eta/2)}/(1+e^-l)^{2+eta} sqrt(cosh(eta l)) dl
def F_t(eta):
    return quad(lambda t: (t*(1-t))**(eta/2)*sqrt(cosh(eta*log(t/(1-t)))), [0, mpf(1)/2, 1])
def F_l(eta):
    return quad(lambda l: mp.e**(-l*(1+eta/2))/(1+mp.e**(-l))**(2+eta)*sqrt(cosh(eta*l)), [-mp.inf, 0, mp.inf])
for e in [mpf('0.2'), mpf('-0.2'), mpf('0.5')]:
    print('eta', e, 'F_t-F_l =', mp.nstr(F_t(e)-F_l(e), 3))
