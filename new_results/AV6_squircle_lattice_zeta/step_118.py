code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum
import mpmath
setdps(60)
s = 3
K = 1200
U2 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        U2 += 4*power(k,-s-2)*(a_+b_)**2
# axes: (|m|=k,n=0),(0,|n|=k): L^2 = k^2 ln^2 k: 4 points
from mpmath import diff as mpdiff
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
U2 += 4*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
print("U_2 direct+axes:", nstr(U2, 20))
# tail: shell ~ 8 sigma_22 ~ 8 k^3/3 ln^2 k /3? use fit
xs = []; ys = []
for k in range(K-30, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        sh += 4*(a_+b_)**2
    xs.append([k**3*log(k)**2, k**3*log(k), k**3, k**2*log(k)**2]); ys.append(sh)
A = mpmath.matrix(31, 4); b = mpmath.matrix(31, 1)
for i in range(31):
    for j in range(4): A[i,j] = xs[i][j]
    b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
tail = mpf(0)
for k in range(K+1, 30*K+1):
    tail += power(k,-s-2)*(c[0,0]*k**3*log(k)**2 + c[1,0]*k**3*log(k) + c[2,0]*k**3 + c[3,0]*k**2*log(k)**2)
print("tail:", nstr(tail, 15))
print("U_2 full:", nstr(U2+tail, 20))
print("formula: 5.22518 ; needed by A1-identity: 5.24938")
'''
exec(code)
