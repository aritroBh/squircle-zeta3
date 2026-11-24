
import mpmath as mp
mp.mp.dps = 100
import json

d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
rt2 = mp.sqrt(2); twopi = 2/mp.pi

def law(n):
    return rt2*mp.cos(mp.pi*(3-2*n)/4)*mp.factorial(n)*twopi**n*n**(-mp.mpf('1.5'))

print("=== ratio f_n / law(n), by residue class mod 4 ===")
for r in range(4):
    ns = [n for n in range(120, 181) if n % 4 == r]
    ratios = [f[n]/law(n) for n in ns]
    print(f"class {r}: n={ns[0]}..{ns[-1]}, ratio at n={ns[-1]}: {mp.nstr(ratios[-1], 20)}")

print()
print("=== ratio at selected n ===")
for n in [40, 80, 120, 160, 180]:
    print(f"n={n}: f_n/law = {mp.nstr(f[n]/law(n), 25)}")

print()
print("=== Richardson / fits per class: ratio = 1 + a/n + b/n^2 ===")
# fit ratio_n = 1 + a/n + b/n^2 + c/n^3 on n in [100,180] per class
import numpy as np
for r in range(4):
    ns = [n for n in range(100, 181) if n % 4 == r]
    X = np.array([[1.0, 1.0/n, 1.0/n**2, 1.0/n**3] for n in ns])
    Y = np.array([float(f[n]/law(n)) for n in ns])
    coef, res, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = np.max(np.abs(X@coef - Y))
    print(f"class {r}: const={coef[0]:.12f}  a={coef[1]:.6f}  b={coef[2]:.4f}  c={coef[3]:.3f}  maxres={resid:.2e}")

print()
print("=== compare f_n vs n! H_n directly (exact H_n) ===")
from fractions import Fraction
from math import factorial
N = 182
KK = N//2+2
a_ = [Fraction(1, factorial(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a_[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*factorial(n-2*k)) for k in range(n//2+1))
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]
for n in [20, 60, 100, 180]:
    v1 = f[n]/(mp.factorial(n)*H[n])
    print(f"n={n}: f_n/(n! H_n) = {mp.nstr(v1, 25)}   (1 - v1 = {mp.nstr(1-v1, 4)})")
print()
print("=== Darboux check on H_n itself: H_n/lawH(n), lawH = law/n! ===")
for n in [40, 80, 120, 180]:
    lawH = rt2*mp.cos(mp.pi*(3-2*n)/4)*twopi**n*n**(-mp.mpf('1.5'))
    print(f"n={n}: H_n/lawH = {mp.nstr(H[n]/lawH, 25)}")
