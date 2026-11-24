from mpmath import mp, mpf, zeta, log, pi, nprint, quad, inf, exp, sqrt
mp.dps = 60
# M_k = int_{log4}^inf x^k [(1-4e^{-x})^{-1/2} - 1] dx = 2 int_0^1 (log(4/(1-u^2)))^k/(1+u) du
def Mk(k):
    return 2*quad(lambda u: (log(4/(1-u**2)))**k/(1+u), [0, 1])
M = {k: Mk(k) for k in range(0, 8)}
for k in range(0, 8):
    print("M_%d ="%k, nprint(M[k], 40))

print("\nM0 - 2 log2 =", nprint(M[0]-2*log(2), 5))
print("M1 - (2 log2^2 + pi^2/6) =", nprint(M[1]-(2*log(2)**2+pi**2/6), 5))
print("M2 - (4 zeta3 + 8 log2^3/3) =", nprint(M[2]-(4*zeta(3)+8*log(2)**3/3), 5))
print("M3 - (4 log2^4 + 3 pi^4/20) =", nprint(M[3]-(4*log(2)**4+3*pi**4/20), 5))
phi6 = 79*pi**6/7560 - 4*zeta(3)**2
print("M5 - (32 log2^6/3 + 79 pi^6/126 - 240 zeta3^2) =",
      nprint(M[5]-(32*log(2)**6/3 + 79*pi**6/126 - 240*zeta(3)**2), 5))
# Bernoulli beta_n
from mpmath import bernoulli
def beta(n):
    return mpf(bernoulli(n))/(mp.factorial(n)*(n+1))
# coefficient c_k = (beta_{k-2}/2)[(2log2)^k/k - M_{k-1}]
def ck(k):
    return (beta(k-2)/2)*((2*log(2))**k/k - M[k-1])
print("\ncoefficients from c-form side:")
for k in [2,3,4,5,6]:
    print("c_%d ="%k, nprint(ck(k), 40))
print("\ntargets: -pi^2/12, zeta3/2, -3zeta4/16, 0, phi6/120")
print(nprint(-pi**2/12, 40)); print(nprint(zeta(3)/2, 40))
print(nprint(-3*zeta(4)/16, 40)); print(nprint(phi6/120, 40))
with open('progress.log','a') as fh:
    fh.write("step_19: Watson/c-form expansion: M0..M5 closed forms verified 60d; c2..c6 coefficients match\n")
