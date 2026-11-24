# step_23.py : comprehensive verification (fast: direct power-series exponentiation)
import math, pickle, time
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta
mp.dps = 70
L2 = log(2); z2 = pi**2/6; z3 = zeta(3); z4 = pi**4/90; z5 = zeta(5); z6 = pi**6/945

K = 30
# log E_g = -L2 z + Σ_{k≥2} c_k z^k, c_k = [(-1)^k(2^k-1)+1] ζ(k)/(k 2^k)
zv = {2: z2, 3: z3, 4: z4, 5: z5, 6: z6}
def zet(k):
    if k in zv: return zv[k]
    return zeta(k)
c = [mpf(0)]*(K+1)
c[1] = -L2
for k in range(2, K+1):
    c[k] = ((-1)**k*(2**k-1)+1)*zet(k)/(k*2**k)
# exp: e = exp(c1 z + c2 z² + ...)
eks = [mpf(0)]*(K+1); eks[0] = mpf(1)
for n in range(1, K+1):
    eks[n] = sum(c[j]*eks[n-j] for j in range(1, n+1))/n
    eks[n] *= 1
# wait: standard recurrence: e_n = (1/n) Σ_{j=1}^n j c_j e_{n-j}
eks = [mpf(0)]*(K+1); eks[0] = mpf(1)
for n in range(1, K+1):
    eks[n] = sum(j*c[j]*eks[n-j] for j in range(1, n+1))/n
mt = [sum(eks[:k+1]) for k in range(K+1)]
print("e_k:", [nstr(eks[k],12) for k in range(9)])
print("m̃_k:", [nstr(mt[k],12) for k in range(9)])
def poly(k, N):
    L = log(mpf(N)); s = L**(k+1)/mpf(k+1)
    for i in range(1, k+2):
        s += mpf(math.factorial(k))/mpf(math.factorial(k+1-i))*mt[i]*L**(k+1-i)
    return s

d = pickle.load(open("E_values_new.pkl","rb"))
print("=== C(n) at O(e^2) ===")
for k,v in d.items():
    n = mpf(k); C = mpf(v[2]); e = 1/n; N = n*sqrt(2)
    p1,p2,p3 = poly(1,N),poly(2,N),poly(3,N)
    Cpred = p1 + e*(p1-p2) - e/(2*sqrt(2)) + e**2*(mpf('0.75')*p1 - mpf('1.25')*p2 + mpf('0.5')*p3 - 3/(2*sqrt(2)))
    print("n=%s: C=%s  diff(C-Cpred)=%s" % (k, nstr(C,36), nstr(C-Cpred,10)))

def R4pred(n):
    N = mpf(n)*sqrt(2); p1,p2,p3 = poly(1,N),poly(2,N),poly(3,N)
    return -3*sqrt(2)/4 - mpf('0.75')*p1 + mpf('1.25')*p2 - mpf('0.5')*p3
b3 = mpf('0.5') - L2/2
c3 = L2**2/4 - L2/2 + 1 + z2
d3 = -L2**3/24 + L2**2/8 - L2/2 - L2*z2/2 + z2/2 - z3/2 + 1 + 3*sqrt(2)/2
def R3(n):
    L = log(mpf(n))
    return L**3/3 + b3*L**2 + c3*L + d3
print()
print("=== residual chain ===")
for k,v in d.items():
    n = mpf(k); E = mpf(v[0]); resid = mpf(v[1]); L = log(n)
    ratio = resid/(L**3/n)
    R3num = n*resid; d3diff = R3num - R3(n); R4 = d3diff*n
    print("n=%s: resid=%s  ratio=%s" % (k, nstr(resid,22), nstr(ratio,10)))
    print("   R3num-R3pred=%s   R4pred=%s   diff=%s   (R5/n est: (log^5)/n^2=%s)" % (nstr(d3diff,14), nstr(R4pred(n),14), nstr(R4-R4pred(n),8), nstr(L**5/n**2,6)))
