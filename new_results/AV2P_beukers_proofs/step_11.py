import sys
sys.path.insert(0, '../AV2_beukers_squircle')
from mpmath import (mp, mpf, nsum, inf, pi, zeta, log, sqrt, quad, gamma, psi,
                    exp, nprint, loggamma, digamma)
mp.dps = 50

def phi(w):
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

# C2 target
C2_target = -pi**2/12
# numerical value of the integral (AV2: -0.82246703342411321823620758332301259461)
C2_num = quad(lambda w: (phi(w)-1)/w**3, [0, 1, inf])
print("C2 integral:", nprint(C2_num, 40))
print("C2 - (-pi^2/12):", nprint(C2_num - C2_target, 5))

# Idea: phi(w) = w^2 2^{-2w} E(w) [w^2 log2 + w/2 + (w^2/2)(psi(w+1/2)-psi(w))], E(w)=B(w,1/2)
def E(w):
    w = mpf(w)
    return exp(loggamma(w) + log(mpf('0.5')) - loggamma(w+mpf('0.5')))
# check decomposition
w = mpf('1.3')
chk = w**2 * 2**(-2*w) * E(w) * (w**2*log(2) + w/2 + (w**2/2)*(digamma(w+mpf('0.5'))-digamma(w)))
print("phi decomposition residual:", nprint(phi(w) - chk, 5))

# So C2 = int [phi(w)-1]/w^3. Split the w^2 log 2 part:
# TermL = int_0^inf 2^{-2w} E(w) log(2)/w dw,  but E(w)~1/w so ~log2/w^2 at 0: divergent; keep combined.
# Full: phi-1 = 2^{-2w}E(w)[w^2 log2 + w/2 + (w^2/2)D2(w)] - 1, D2 = psi(w+1/2)-psi(w)
# Group as: [2^{-2w}E(w) w/2 - 1] + 2^{-2w}E(w) w^2 log2 + 2^{-2w}E(w)(w^2/2)D2(w)
# each /w^3:
#  A = [2^{-2w}E(w) w/2 - 1]/w^3
#  B = 2^{-2w}E(w) log2 / w
#  C = 2^{-2w}E(w) D2(w)/(2w)
# check A,B,C integrability: A at 0: 2^{-2w}E(w)w/2 -> 1, expansion?
# E(w) = sqrt(pi) G(w)/G(w+1/2); series at 0:
import mpmath as m
m.mp.dps = 30
s = lambda f: m.taylor(f, 0, 4)
Es = s(lambda w: E(w))
print("E(w) taylor:", Es)
