code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
setdps(60)
def zetap(a): return mpdiff(lambda x: zeta(x), a)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
s = 4
c1 = mpf('1.245951431164482619524932648613589408233')
# Z''(1) = sum T [s^2 (lnS - w)^2 - 2s ln S + 2s w - s(M/S - w^2)]
K = 1500
Zpp = mpf(0)
for k in range(1, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        lS = log(k)
        sh += 4*(s**2*(lS-w)**2 - 2*s*lS + 2*s*w - s*(M/k - w*w))
    # axes: |m|=k,n=0: w = ln k, M/S = ln^2 k, (lnS-w)=0: contributes 2*(-2s ln k + 2s ln k - s(ln^2k - ln^2k)) = 0
    Zpp += power(k, -s)*sh
# tail ~ k^{1-s} ln^2 k : use EM estimate: shell average ~ 4k * typical:
# typical bracket ~ s^2 (ln k)^2 * (variance of w) ... compute tail via last-shell trend
# estimate tail by fitting shell(k) ~ k*(a ln^2 k + b ln k + c):
sh3 = [ (kk, None) for kk in range(K-200, K+1)]
# just do EM with the integral model: bracket_avg ~ s^2 ln^2 k * Var(w) with Var(w) ~ (ln k)^2/... 
# Simpler: fit shell(k) = k*(a ln^2 k + b ln k + c) from last 30 shells
import mpmath
xs = []
ys = []
for k in range(K-30, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        lS = log(k)
        sh += 4*(s**2*(lS-w)**2 - 2*s*lS + 2*s*w - s*(M/k - w*w))
    xs.append([k*log(k)**2, k*log(k), k])
    ys.append(sh)
A = mpmath.matrix(31, 3); b = mpmath.matrix(31, 1)
for i in range(31):
    for j in range(3): A[i,j] = xs[i][j]
    b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
tail = mpf(0)
for k in range(K+1, 30*K+1):
    model = c[0,0]*k*log(k)**2 + c[1,0]*k*log(k) + c[2,0]*k
    tail += power(k, -s)*model
# beyond 30K: model*k^{-s} ~ ln^2 k/k^{s-1} : for s=4 ~ ln^2/k^2: add integral estimate
tail += c[0,0]*(zetapp(s-1) - fsum([power(m,1-s)*log(m)**2 for m in range(1,30*K+1)])) 
Zpp += tail
print("Z''(1) =", nstr(Zpp, 25))
print("c_2 = c_1 + Z''/2 =", nstr(c1 + Zpp/2, 25))
print("fit: 1.573849628813923486848")
'''
exec(code)
