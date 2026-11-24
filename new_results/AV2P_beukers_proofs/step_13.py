import sys
sys.path.insert(0, '../AV2_beukers_squircle')
from mpmath import (mp, mpf, nsum, inf, pi, zeta, log, sqrt, quad, gamma, psi,
                    exp, nprint, loggamma, digamma)
mp.dps = 60

def phi(w):
    w = mpf(w)
    if abs(w) < mpf('0.05'):
        # Taylor: 1 - zeta(3) w^3 + pi^4/40 w^4 + phi5 w^5 + phi6 w^6 ...
        z3 = zeta(3); p4 = pi**4/40
        p5 = 3*zeta(2)*zeta(3)-9*zeta(5)
        p6 = 79*pi**6/7560 - 4*zeta(3)**2
        return 1 - z3*w**3 + p4*w**4 + p5*w**5 + p6*w**6
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

def E(w):
    w = mpf(w)
    return exp(loggamma(w) + log(sqrt(pi)) - loggamma(w+mpf('0.5')))

# check decomposition phi(w) = 2^{-2w} E(w) [(w^2/2)(psi(w+1/2)-psi(w)) + w^2 log2 + w/2]
for w0 in ['1.3','0.37','2.9','0.11']:
    w = mpf(w0)
    chk = 2**(-2*w)*E(w)*((w**2/2)*(digamma(w+mpf('0.5'))-digamma(w)) + w**2*log(2) + w/2)
    print(f"phi decomp residual at {w0}:", nprint(phi(w)-chk, 4))

# C2 integral, careful split
C2_num = quad(lambda w: (phi(w)-1)/w**3, [0, 1, 10, 100, inf])
print("C2 integral:", nprint(C2_num, 45))
print("C2 - (-pi^2/12):", nprint(C2_num + pi**2/12, 5))
