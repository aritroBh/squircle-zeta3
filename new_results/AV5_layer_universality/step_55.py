import mpmath
from mpmath import mp, mpf, sqrt, zeta, log, cos, sin, gamma, euler, nstr, factorial
mp.dps = 100
def psi_m(m, n):
    if m==0: return -euler + sum(mpf(1)/j for j in range(1,n))
    return (-1)**(m+1)*factorial(m)*(zeta(m+1) - sum(mpf(1)/j**(m+1) for j in range(1,n)))
Aser_log = []
for k in range(1,11):
    ak = (2*(-1)**k*psi_m(k-1,2) - (-2)**k*psi_m(k-1,3))/factorial(k)
    Aser_log.append(ak)
def exp_series(coeffs, N):
    b = [mpf(1)]
    for n in range(1,N+1):
        bn = sum(k*coeffs[k-1]*b[n-k] for k in range(1,n+1))/n
        b.append(bn)
    return b
b = exp_series(Aser_log, 10)
print("log(A/4) = sum a_k tau^k:")
for k,a in enumerate(Aser_log,1): print(f"  a{k} = {nstr(a,38)}")
print("A(tau)/4 coefficients:")
for k,bk in enumerate(b): print(f"  [tau^{k}] {nstr(bk,38)}")
tau = mpf('0.01')
Adirect = 4*gamma(2-tau)**2/gamma(3-2*tau)
Aser = 4*sum(b[k]*tau**k for k in range(11))
print("direct:", nstr(Adirect,55)); print("series:", nstr(Aser,55)); print("diff:", nstr(Adirect-Aser,5))
