# P2N4: the explicit constants of the proof, in exact rational arithmetic (every rounding is upward).
from fractions import Fraction as Fr
import math
beta = Fr(9, 8)                     # 1/(1-x) <= 9/8 for x = 1/n <= 1/9
L0 = Fr(21972, 10000)               # log 9 = 2.19722... >= L0
L0up = Fr(21973, 10000)             # log 9 <= L0up
assert L0 <= math.log(9) <= L0up
A = Fr(215, 100)                    # a2 < zeta(2)+1/2 < 2.15
assert math.pi**2/6 + 0.5 < 2.15
a3max = Fr(31, 10)*5**3             # Cauchy at 0, radius 1/5: |a3| <= 3.1*5^3
b4max = Fr(31, 10)*5**4             # |b4| <= 3.1*5^4
K5 = Fr(31, 10)*10**5               # |rho_5| <= 3.1*10^5 eps^5 (Lagrange + Cauchy radius 1/10)
s1 = (L0up+1)/9                     # x*lam <= (log 9 + 1)/9
s2 = (L0up+1)**2/9                  # x*lam^2 <= (log 9 + 1)^2/9
kappa = 1 + 1/L0                    # lam <= L+1 <= kappa L
b5, b6 = beta**5, beta**6
# coefficients of x^5 * (c0 + c1 lam + c2 lam^2) * z, piece by piece
pieces = {
 'P1 (Taylor of (1-x)^-k)': (2*A*4*b5, 2*A*6*b5, 2*A*2*b5),
 'P1 (theta*eps^5 lam^3 - L^3 x^5)': (0, 0, A/3*(5*b6*s1 + 3 + b6*s2)),
 'P2': (2*a3max*6*b5, 2*a3max*4*b5, 2*a3max*b5/2),
 'P3': (2*b4max*4*b5, 2*b4max*b5, 0),
 'P4': (2*K5*b5, 0, 0),
 'P5': (16*Fr(9, 7), 0, 0),
 'P6 m=2': (2*Fr(9, 8)/Fr(5, 6), 0, 0),
 'P6 m=3': (2*Fr(4, 3)/(3-Fr(4, 9)), 0, 0),
 'P6 m>=4': (2*Fr(9, 31)*beta, 0, 0),
}
C0 = C1 = C2 = 0
for k, (c0, c1, c2) in pieces.items():
    print("%-36s c0=%12.4f c1=%10.4f c2=%9.4f" % (k, float(c0), float(c1), float(c2)))
    C0 += c0; C1 += c1; C2 += c2
print("totals: c0=%.4f c1=%.4f c2=%.4f" % (float(C0), float(C1), float(C2)))
# convert to L^2: 1 <= L^2/L0^2, lam <= kappa L <= kappa L^2/L0, lam^2 <= kappa^2 L^2
Cp = C0/L0**2 + C1*kappa/L0 + C2*kappa**2
print("C' (coefficient of z L^2 x^5) =", float(Cp))
C4 = A/3 + Cp/L0
print("C4 (coefficient of z L^3 x^5) =", float(C4))
print("used in result.md: C' <= 2.5e5, C4 <= 1.2e5:", Cp <= 250000, C4 <= 120000)
print("aux: beta^5=%.5f beta^6=%.5f s1=%.5f s2=%.5f kappa=%.5f" % tuple(map(float, (b5, b6, s1, s2, kappa))))
