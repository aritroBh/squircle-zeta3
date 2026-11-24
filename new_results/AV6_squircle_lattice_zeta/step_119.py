code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff, pi
import mpmath
setdps(60)
s = 3
def h(x):
    return -x*log(x) - (1-x)*log(1-x)
Ch = 4*(mpf(5)/6 - pi**2/18)   # 4 * integral h^2 = 4*(5/6 - pi^2/18)
K = 1200
# exact shells for k in [K-100, K]
xs = []; ys = []
for k in range(K-100, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        w = (a_+b_)/k
        sh += 4*(log(k)-w)**2
    xs.append([k, log(k), mpf(1), 1/k]); ys.append(sh)
A = mpmath.matrix(len(xs), 4); b = mpmath.matrix(len(xs), 1)
for i in range(len(xs)):
    for j in range(4): A[i,j] = xs[i][j]
    b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
print("shell model: a=", nstr(c[0,0], 15), " (4C_h =", nstr(Ch, 15), ")")
print("             b=", nstr(c[1,0], 12), " c=", nstr(c[2,0], 12), " d=", nstr(c[3,0], 12))
def HT(sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
def HTln(sig):
    return -mpdiff(lambda x: zeta(x), sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
tail = c[0,0]*HT(s-1) + c[1,0]*HTln(s) + c[2,0]*HT(s) + c[3,0]*HT(s+1)
print("tail:", nstr(tail, 15))
# A1 partial
A1 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        w = (a_+b_)/k
        A1 += 4*power(k,-s)*(log(k)-w)**2
print("A1 full:", nstr(A1+tail, 20))
'''
exec(code)
