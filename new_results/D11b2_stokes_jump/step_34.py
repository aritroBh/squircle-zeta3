# ======== STEP 2d: fits — exponent, constant, divergence ========
from mpmath import mp, mpf, pi, factorial, log
import json
mp.dps = 120
F = [mpf(v) for v in json.load(open('step2_num200.json'))['f200']]
def pat(n): return [-1, 1, 1, -1][n % 4]
R = {n: F[n]/(factorial(n)*(2/pi)**n) for n in range(2, 65)}
# 1) effective exponent nu_n from ratio of successive same-parity R
print("effective exponent nu (from R_n/R_{n+2}): expect -> 3/2")
for n in [30, 40, 50, 58, 59]:
    nu = log(abs(R[n]/R[n+2]))/log((n+2)/n)
    print(f"  n={n:3d}  nu = {mp.nstr(nu, 12)}")
# 2) constant via Richardson on C_n = R_n n^{3/2} / pat(n) = 1 + a/n + b/n^2 + c/n^3
C = {n: R[n]*n**mpf('1.5')/pat(n) for n in range(10, 64)}
def rich(vals, ns, p):  # fit C = K + sum a_i n^-i, i=1..p, return K
    import numpy as np
    M = []
    for n in ns: M.append([n**-i for i in range(1, p+1)] + [1.0])
    A = mp.matrix([[mpf(M[r][c]) for c in range(p+1)] for r in range(p+1)])
    b = mp.matrix([vals[n] for n in ns])
    x = (A**-1)*b
    return x[p]
print("\nconstant extrapolation (fit 1 + a/n + b/n^2 (+c/n^3)):")
for ns in [[56, 58, 60], [54, 58, 62], [57, 59, 61], [53, 57, 61], [50, 54, 58, 62], [51, 55, 59, 63]]:
    print(f"  ns={ns}: K = {mp.nstr(rich(C, ns, len(ns)-1), 15)}")
# 3) divergence: |f_n|^{1/n}
print("\n|f_n|^{1/n} vs (2/pi)(n/e) =", mp.nstr(2/pi/mp.e, 8), "* n:")
for n in [20, 40, 60]:
    print(f"  n={n}: |f_n|^(1/n) = {mp.nstr(abs(F[n])**(1/n), 10)}   (2/pi)(n/e) = {mp.nstr(2*n/pi/mp.e, 10)}")
