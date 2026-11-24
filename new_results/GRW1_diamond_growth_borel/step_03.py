
import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial

# Exact rational H_n for square: H(z) = 2 e^{-z/2} sqrt(cosh z) = sqrt(2)*sqrt(1+e^{-2z})
N = 500
K = N//2+2
a = [Fraction(1, factorial(2*k)) for k in range(K)]
s = [Fraction(0)]*K; s[0] = Fraction(1)
for k in range(1, K):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
H = [Fraction(0)]*(N+1)
for n in range(N+1):
    H[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*factorial(n-2*k)) for k in range(n//2+1))

# check H(z) = sqrt(2) sqrt(1+e^{-2z}) coefficient identity: H_n should equal [z^n] sqrt(2)*sqrt(1+e^{-2z})
# expand: sqrt(2) sum_k binom(1/2,k) e^{-2kz}; symbolic check on the SERIES level:
# sqrt(1+w) = sum c_k w^k ; H(z) = sqrt(2) sum_k c_k e^{-2kz}  => for Re z>0.
# Compare H(z) value vs sqrt(2)*sqrt(1+e^{-2z}) at some z:
def Hval(z):
    return 2*mp.exp(-z/2)*mp.sqrt(mp.cosh(z))
for z in [mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('-0.2')]:
    v1 = Hval(z); v2 = mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
    print(f"z={z}: H={mp.nstr(v1,25)} sqrt2*sqrt(1+e^-2z)={mp.nstr(v2,25)} diff={mp.nstr(abs(v1-v2),3)}")

def Hb0(x):
    tot = mp.mpf(0)
    for m in range(N+1):
        tot += mp.mpf(H[m].numerator)/mp.mpf(H[m].denominator) * x**m/mp.factorial(m)
    return tot

# S(x) := sqrt(2) sum_k c_k J_0(2 sqrt(2 k x)), c_k = binom(1/2,k)
ck = [mp.mpf(1)]
for k in range(1, 20001):
    ck.append(ck[-1]*(mp.mpf('0.5')-(k-1))/k)
def S(x, kmax=20000):
    tot = mp.mpf(0)
    for k in range(kmax+1):
        tot += ck[k]*mp.besselj(0, 2*mp.sqrt(2*k*x))
    return mp.sqrt(2)*tot

print()
for x in [5, 20, 80]:
    v1 = Hb0(mp.mpf(x))
    v2 = S(mp.mpf(x), 20000)
    print(f"x={x}: H^[0](x)={mp.nstr(v1,25)}")
    print(f"        S(x)   ={mp.nstr(v2,25)}  diff={mp.nstr(abs(v1-v2),3)}  (sqrt2={mp.nstr(mp.sqrt(2),12)})")
