"""SRF4 coefficient engine (independent of SRF1/SRF3 engines).
Gamma-variable identity:  Gamma(3-2tau) T(tau) = E[(g1g2g3)^(-tau) (sum g_i^(2tau))^(1/2)],  g_i iid Exp(1),
expanded exactly: (sum g^{2tau})^{1/2} = sqrt3 sum_n binom(1/2,n) 3^-n (sum X_i)^n, X_i = g_i^{2tau}-1,
E[g^-tau X^m] = sum_k [x^k](e^{-x}(e^{2x}-1)^m) * Gamma^{(k)}(1) tau^k.
Outputs a_k=[tau^k]S/S(1), b_n=[d^n]S/S(1).  Usage: python coeffs.py K dps"""
import sys, json
from mpmath import mp, mpf, zeta, euler, binomial, factorial, sqrt, pi, log, nstr
K = int(sys.argv[1]); mp.dps = int(sys.argv[2])
def mul(a, b):
    c = [mpf(0)]*(K+1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(K+1-i): c[i+j] += x*b[j]
    return c
def exp_ser(a):  # exp of series with a[0]=0
    e = [mpf(0)]*(K+1); e[0] = mpf(1)
    for n in range(1, K+1):
        e[n] = sum(k*a[k]*e[n-k] for k in range(1, n+1))/n
    return e
# gamma_k = [x^k] Gamma(1+x)
lg = [mpf(0)]*(K+1); lg[1] = -euler
for k in range(2, K+1): lg[k] = (-1)**k*zeta(k)/k
gam = exp_ser(lg)
EL = [factorial(k)*gam[k] for k in range(K+1)]           # E[(log g)^k]
# phi_m series
em1 = [mpf(0)] + [mpf(2)**k/factorial(k) for k in range(1, K+1)]   # e^{2x}-1
emx = [(-1)**k/factorial(k) for k in range(K+1)]                   # e^{-x}
A = []; P = [mpf(1)] + [mpf(0)]*K
for m in range(K+1):
    phi = mul(emx, P)
    A.append([phi[k]*EL[k] for k in range(K+1)])
    P = mul(P, em1)
# B(y,tau) = sum_m A_m y^m/m!  as list over m of tau-series
B = [[A[m][k]/factorial(m) for k in range(K+1)] for m in range(K+1)]
def bimul(X, Y):
    Z = [[mpf(0)]*(K+1) for _ in range(K+1)]
    for m1 in range(K+1):
        for m2 in range(K+1-m1):
            # tau-degree >= m1+m2
            x, y, z = X[m1], Y[m2], Z[m1+m2]
            for k1 in range(m1, K+1):
                xv = x[k1]
                if xv == 0: continue
                for k2 in range(m2, K+1-k1): z[k1+k2] += xv*y[k2]
    return Z
B3 = bimul(bimul(B, B), B)
F = [mpf(0)]*(K+1)
for n in range(K+1):
    c = sqrt(3)*binomial(mpf(1)/2, n)/mpf(3)**n*factorial(n)
    for k in range(K+1): F[k] += c*B3[n][k]
# 1/Gamma(3-2tau) = 1/((2-2t)(1-2t)Gamma(1-2t))
G1m2t = [gam[k]*(-2)**k for k in range(K+1)]
G3 = mul(mul([mpf(2), mpf(-2)] + [mpf(0)]*(K-1), [mpf(1), mpf(-2)] + [mpf(0)]*(K-1)), G1m2t)
inv = [mpf(0)]*(K+1); inv[0] = 1/G3[0]
for n in range(1, K+1): inv[n] = -sum(G3[k]*inv[n-k] for k in range(1, n+1))/G3[0]
T = mul(F, inv)
a = mul(mul([mpf(1), mpf(-1)] + [mpf(0)]*(K-1), [mpf(1), mpf(-1)] + [mpf(0)]*(K-1)), T)
a = [2/sqrt(3)*x for x in a]
# b_n = sum_k a_k [d^n](d/(1+d))^k = sum_k a_k (-1)^(n-k) C(n-1,k-1)
b = [a[0]] + [sum(a[k]*(-1)**(n-k)*binomial(n-1, k-1) for k in range(1, n+1)) for n in range(1, K+1)]
json.dump({"K": K, "dps": mp.dps, "T": [str(x) for x in T], "a_tau": [str(x) for x in a], "b_d": [str(x) for x in b]},
          open(f"coeffs_K{K}_dps{mp.dps}.json", "w"))
print("T0 =", nstr(T[0], 20), " sqrt3/2 =", nstr(sqrt(3)/2, 20))
z3 = zeta(3)
ref = [1, 1, 1-pi**2/9, 1+pi**2/9-16*z3/9, 1-2*pi**2/9+32*z3/9-16*pi**4/405]
print("d-series vs paper thm:oct:", [nstr(b[n]-ref[n], 3) for n in range(5)])
