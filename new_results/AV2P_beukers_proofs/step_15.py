import sys
sys.path.insert(0, '../AV2_beukers_squircle')
from mpmath import (mp, mpf, inf, pi, zeta, log, sqrt, quad, exp, nprint,
                    loggamma, digamma)
mp.dps = 80

def phi(w):
    w = mpf(w)
    if abs(w) < mpf('0.05'):
        z3 = zeta(3); p4 = pi**4/40
        p5 = 3*zeta(2)*zeta(3)-9*zeta(5)
        p6 = 79*pi**6/7560 - 4*zeta(3)**2
        p7 = 0  # phi_7? need for w<0.05: w^7 ~ 8e-10, * coeff maybe big; compute exactly below
        return 1 - z3*w**3 + p4*w**4 + p5*w**5 + p6*w**6
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

# check Taylor accuracy at w=0.049
w = mpf('0.049')
direct = w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))
print("taylor vs direct at 0.049:", nprint(phi(w)-direct, 4))

# integrate in pieces
parts = [(0, mpf('0.5')), (mpf('0.5'), 2), (2, 8), (8, 30), (30, 100), (100, inf)]
tot = mpf(0)
for a,b in parts:
    v = quad(lambda w: (phi(w)-1)/w**3, [a, b])
    tot += v
print("C2 =", nprint(tot, 70))
print("C2 + pi^2/12 =", nprint(tot + pi**2/12, 5))
