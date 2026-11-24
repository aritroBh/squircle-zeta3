# Planner check of SRF1: Borel rate of the octahedron-end tau-series from its stored 60 coefficients (Hankel ratios), vs 2/sqrt(pi^2+log^2 2).
import json
from mpmath import mp, mpf, factorial, sqrt, log, pi, nstr, atan2
mp.dps = 120
a = [mpf(x) for x in json.load(open('SRF1_n3_K60.json'))['tau']]
b = [a[k]/factorial(k) for k in range(len(a))]
u = [b[k+1]**2 - b[k]*b[k+2] for k in range(len(b)-2)]
est = [sqrt(abs(u[k+1]/u[k])) for k in range(len(u)-1)]
# remove the k^(2beta) drift: A_k = est_k ~ A(1 + beta'/k): Richardson on k*est
ks = list(range(len(est)))
rich = [(k+1)*est[k+1] - k*est[k] for k in ks[30:-1]]
print("Hankel-ratio estimates k=40..56:", [nstr(est[k], 6) for k in range(40, len(est), 4)])
print("Richardson (1st order) last 4:", [nstr(r, 6) for r in rich[-4:]])
print("predicted A = 2/sqrt(pi^2+log^2 2) =", nstr(2/sqrt(pi**2 + log(2)**2), 8), "; theta =", nstr(atan2(pi, -log(2)), 8))
# Cube end: residue of V at eps = -1. Planner identification: R = I - 2, I = int_0^1 (sqrt(1+a^4)-1)/a^2 da = 1 - 2 pi^(3/2)/Gamma(1/4)^2.
from mpmath import quad, gamma
I = quad(lambda t: (sqrt(1 + t**4) - 1)/t**2, [0, 1])
R = -1 - 2*pi**mpf(1.5)/gamma(mpf(1)/4)**2
V = json.load(open('EV_coeffs_K30.json'))['V']
print("I - (1 - 2pi^1.5/Gamma(1/4)^2) =", nstr(I - (1 - 2*pi**mpf(1.5)/gamma(mpf(1)/4)**2), 5))
print("R =", nstr(R, 15), "; (-1)^k V_k, k = 26..30:", [nstr(mpf(V[k])*(-1)**k, 10) for k in range(26, 31)])
