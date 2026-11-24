"""v_taylor_exact.py -- verifier's EXACT check (Fractions) that the HPL forms B1, B2 (result.md, typed by hand)
have Taylor coefficients equal to the harmonic-number formulas of Theorem A (typed by hand from result.md):
  mu_{1,n} = D_n/(n(2n-1)) - 1/(2n^2),
  mu_{2,n} = z2/(4n^2-1) - D_n^2/(2n(2n-1)) + D_n/(2n^2) + H2_{n-1}/(2n(4n^2-1)) - H2_{2n-1}/(2n(2n-1))
             - (4n^2-2n-1)/(4n^3(2n+1)),   D_n = H_{2n-1} - H_{n-1}.
B1 (HPL): M1 = -(1-c)[H_{0,1} + H_{1,-1}] + (1+c)[H_{0,-1} - H_{-1,1}]
B2 (HPL): M2 = z2 Phi(c) - Li2(c^2)/2 + S(c) + S(-c)   [z2/2 - z2 q artanh c = z2 Phi]
   S(c) = (1-c)[H001 + H01m1 - H1m11] + q H101 - ((1-c)^2/(2c)) H10m1,  q = (1-c^2)/(2c).
Own HPL power-series code (exact rationals); S(-c) is evaluated by substituting -c in the series directly.
Also checks the M0 formula mu_{0,n} = -1/(n(2n-1)) against -(1+c)log(1+c)-(1-c)log(1-c).
"""
import sys
from fractions import Fraction as Fr

NN = int(sys.argv[1]) if len(sys.argv) > 1 else 250      # check n <= NN  (powers c^(2n))
NT = 2*NN + 4


def ser_H(word):
    a = [Fr(0)]*(NT + 1)
    last = word[-1]
    for k in range(1, NT + 1):
        a[k] = Fr(1, k) if last == 1 else Fr((-1)**(k + 1), k)
    for letter in reversed(word[:-1]):
        b = [Fr(0)]*(NT + 1)
        if letter == 0:
            for k in range(1, NT + 1):
                b[k] = a[k]/k
        else:
            acc = Fr(0)
            for m in range(1, NT):
                acc = (-acc + a[m]) if letter == -1 else (acc + a[m])
                b[m + 1] = acc/(m + 1)
        a = b
    return a


def neg(s):                      # series of f(-c)
    return [v*(-1)**k for k, v in enumerate(s)]


def mulpoly(s, poly):            # poly: dict exponent -> Fraction (exponents may be negative)
    out = {}
    for e, co in poly.items():
        for k, v in enumerate(s):
            if v:
                out[k + e] = out.get(k + e, Fr(0)) + co*v
    return out


def addto(d, src, f=Fr(1)):
    for k, v in src.items():
        d[k] = d.get(k, Fr(0)) + f*v


H = {w: ser_H(w) for w in [(0, 1), (1, -1), (0, -1), (-1, 1), (0, 0, 1), (0, 1, -1), (1, -1, 1), (1, 0, 1), (1, 0, -1)]}
Hn = {w: neg(s) for w, s in H.items()}
# B1
M1 = {}
addto(M1, mulpoly([x + y for x, y in zip(H[(0, 1)], H[(1, -1)])], {0: Fr(-1), 1: Fr(1)}))
addto(M1, mulpoly([x - y for x, y in zip(H[(0, -1)], H[(-1, 1)])], {0: Fr(1), 1: Fr(1)}))


def S_series(Hs, sign):          # S(x) at x = sign*c, expressed as a series in c
    x1 = {0: Fr(1), 1: Fr(-sign)}               # 1 - x
    qx = {-1: Fr(sign, 2), 1: Fr(-sign, 2)}       # (1 - x^2)/(2x)
    sq = {-1: Fr(-sign, 2), 0: Fr(1), 1: Fr(-sign, 2)}   # -(1-x)^2/(2x) = -(1/(2x) - 1 + x/2)
    out = {}
    addto(out, mulpoly([p + q - r for p, q, r in zip(Hs[(0, 0, 1)], Hs[(0, 1, -1)], Hs[(1, -1, 1)])], x1))
    addto(out, mulpoly(Hs[(1, 0, 1)], qx))
    addto(out, mulpoly(Hs[(1, 0, -1)], sq))
    return out


R2 = {}
addto(R2, S_series(H, 1)); addto(R2, S_series(Hn, -1))
# -Li2(c^2)/2 : Li2(c^2) = sum c^(2k)/k^2
for k in range(1, NN + 2):
    R2[2*k] = R2.get(2*k, Fr(0)) - Fr(1, 2*k*k)


def Hh(r, m):
    return sum((Fr(1, t**r) for t in range(1, m + 1)), Fr(0))


bad = []
negpow = [k for k, v in list(M1.items()) + list(R2.items()) if k <= 0 and v != 0]
odd = [k for k, v in list(M1.items()) + list(R2.items()) if 0 < k <= 2*NN and k % 2 == 1 and v != 0]
h1 = {};
for n in range(1, NN + 1):
    Dn = Hh(1, 2*n - 1) - Hh(1, n - 1)
    mu1 = Dn/(n*(2*n - 1)) - Fr(1, 2*n*n)
    r2 = -Dn**2/(2*n*(2*n - 1)) + Dn/(2*n*n) + Hh(2, n - 1)/(2*n*(4*n*n - 1)) - Hh(2, 2*n - 1)/(2*n*(2*n - 1)) \
        - Fr(4*n*n - 2*n - 1, 4*n**3*(2*n + 1))
    if M1.get(2*n, Fr(0)) != mu1:
        bad.append(('M1', n))
    if R2.get(2*n, Fr(0)) != r2:
        bad.append(('R2', n))
print(f'nonzero coefficients at c^k, k<=0: {negpow}; odd powers up to c^{2*NN}: {len(odd)}')
print(f'B1 (HPL) Taylor coeffs == mu_(1,n) and rational part of B2 (HPL) == r_n for all n <= {NN}:', not bad, bad[:5])
# M0: -(1+c)log(1+c)-(1-c)log(1-c) coefficients vs -1/(n(2n-1))
L1 = ser_H((1,)); Lm1 = ser_H((-1,))   # -log(1-c), log(1+c)
M0 = {}
addto(M0, mulpoly(L1, {0: Fr(1), 1: Fr(-1)}))      # (1-c)(-log(1-c))
addto(M0, mulpoly(Lm1, {0: Fr(-1), 1: Fr(-1)}))    # -(1+c)log(1+c)
print('M0 coeffs == -1/(n(2n-1)) for n <=', NN, ':', all(M0.get(2*n, Fr(0)) == Fr(-1, n*(2*n - 1)) for n in range(1, NN + 1)))
