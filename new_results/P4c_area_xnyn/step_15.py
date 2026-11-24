# step_13.py : verify Tauberian lemma for T_k and the Delta2 series, term by term
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, gamma, inf, quad
import sympy as sp
mp.dps = 60
L2 = log(2); z2 = zeta(2)

def gW(w):
    w = mpf(w); return (1 + w - sqrt(1+w*w))/w

zz = sp.symbols('zz')
Eg_series = sp.series(sp.gamma((1+zz)/2)*sp.gamma(1-zz/2)/sp.sqrt(sp.pi), zz, 0, 8).removeO().expand()
eks = [mpf(float(sp.N(Eg_series.coeff(zz, k), 30))) for k in range(8)]
m = [sum(eks[:k+1]) for k in range(8)]
print("e_k =", [nstr(x,14) for x in eks])
print("m_k =", [nstr(x,14) for x in m])
# closed form checks
print("e1 + log2 =", nstr(eks[1]+L2, 8), " e2 - (L2^2/2+z2/2) =", nstr(eks[2]-(L2**2/2+z2/2),8))

import math
def T_k_num(N, k):
    N = mpf(N)
    return quad(lambda w: log(N/w)**k * gW(w), [0, 1, N])
def T_k_pred(N, k):
    L = log(mpf(N))
    s = L**(k+1)/mpf(k+1)
    for i in range(k+1):
        s += mpf(math.factorial(k))/mpf(math.factorial(k-i)) * m[i] * L**(k-i)
    return s
for N in [mpf(100)*sqrt(2), mpf(10000)*sqrt(2)]:
    for k in [0,1,2,3,4]:
        tn = T_k_num(N,k); tp = T_k_pred(N,k)
        print("N=%s k=%d: T=%s pred=%s diff=%s" % (nstr(N,8), k, nstr(tn,18), nstr(tp,18), nstr(tn-tp,6)))

# Delta2 series vs direct
def D2_direct(n):
    n = mpf(n); e = 1/n; N = n*sqrt(2)
    return quad(lambda w: ((w/N)**e - 1)*log(N/w)*gW(w), [0, 1, N])
def D2_series(n, jmax=8):
    n = mpf(n); e = 1/n; N = n*sqrt(2)
    s = mpf(0)
    for j in range(1, jmax+1):
        s += (-1)**j * e**j/math.factorial(j) * T_k_pred(N, j+1)
    return s
for n in [100, 1000]:
    d2 = D2_direct(n); d2s = D2_series(n)
    print("n=%d: Delta2 direct=%s  series=%s  diff=%s" % (n, nstr(d2,20), nstr(d2s,20), nstr(d2-d2s,8)))
    # also via C' - C0
