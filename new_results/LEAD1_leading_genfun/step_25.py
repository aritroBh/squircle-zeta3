import json, mpmath
mpmath.mp.dps = 410
from mpmath import mpf, log, pi, zeta, gamma, factorial
cw = {int(w): mpf(v) for w, v in json.load(open('cw_dps420.json')).items()}
lam = log(2)
# E(y) = 2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2) = exp(sum s_n y^n/n), s_n = (1-2^{1-n})zeta(n) odd n, zeta(n) even n
# build e_w = [y^w] E(y) exactly in Gamma basis via series: log E = y log2 + sum_{n>=2} s_n y^n/n
NN = 12
import sympy as sp
y = sp.symbols('y')
logE = lam_s = sp.symbols('lam')  # keep symbolic in zetas? simpler: numeric series at 410 dps
# numeric approach with mpmath series via exp of power series
def eta(n): return (1 - mpf(2)**(1-n)) * zeta(n)
s = {n: (zeta(n) if n%2==0 else eta(n)) for n in range(2, NN+1)}
# series coefficients of logE: L_1 = lam, L_n = s_n/n
L = [mpf(0)]*(NN+1)
L[1] = lam
for n in range(2, NN+1): L[n] = s[n]/n
# exp series
E = [mpf(0)]*(NN+1); E[0] = mpf(1)
for n in range(1, NN+1):
    E[n] = sum(k*L[k]*E[n-k] for k in range(1, n+1))/n
# zigzag numbers A_n: sec x + tan x
A = [1,1,1,2,5,16,61,272,1385,7936,50521]
def eps(m):
    return [1,-1,-1,1][m%4]
print("w | c(w) | law | diff")
for w in range(2, 10):
    m = w-2
    law = eps(m)*A[m+1]/(2*factorial(m)) * E[w]
    print(w, mpmath.nstr(cw[w], 25), mpmath.nstr(law, 25), mpmath.nstr(cw[w]-law, 8))
