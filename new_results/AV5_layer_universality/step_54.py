import mpmath, json
from mpmath import mp, mpf, sqrt, zeta, log, cos, sin, gamma, psi, nstr, matrix, lu_solve, factorial
mp.dps = 100
# ---- AREA at diamond end: A = 4 G(2-t)^2/G(3-2t), t = 1-1/p ----
# exact log-coefficients via polygamma at integers
def psi_m(m, n):  # psi^{(m)}(n), integer n>=1
    return (-1)**(m+1)*factorial(m)*(zeta(m+1) - sum(mpf(1)/j**(m+1) for j in range(1,n)))
Aser_log = []
for k in range(1,11):
    ak = (2*(-1)**k*psi_m(k-1,2) - (-2)**k*psi_m(k-1,3))/factorial(k)
    Aser_log.append(ak)
print("log(A/4) coefficients a_k (tau^k):")
for k,a in enumerate(Aser_log,1):
    print(f"  a{k} = {nstr(a,40)}")
# exponentiate
import sympy as sp
t = sp.symbols('t')
logA = sum(sp.Rational(str(float(0)))*t for _ in range(0))  # placeholder
# do exponentiation numerically instead
def exp_series(coeffs, N):
    # coeffs: a1..aN of log; return b0..bN of exp
    b = [mpf(1)]
    for n in range(1,N+1):
        bn = sum(k*coeffs[k-1]*b[n-k] for k in range(1,n+1))/n
        b.append(bn)
    return b
b = exp_series(Aser_log, 10)
print("A(t)/4 = exp(log): coefficients of tau^k:")
for k,bk in enumerate(b):
    print(f"  [tau^{k}] {nstr(bk,40)}")
# verify against direct gamma at t=0.01
tau = mpf('0.01')
Adirect = 4*gamma(2-tau)**2/gamma(3-2*tau)
Aser = 4*sum(b[k]*tau**k for k in range(11))
print("direct:", nstr(Adirect,60)); print("series:", nstr(Aser,60)); print("diff:", nstr(Adirect-Aser,5))
