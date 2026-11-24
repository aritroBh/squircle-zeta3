code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum
import mpmath
setdps(60)
s = 3
K = 1200
A1 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        w = (a_+b_)/k
        A1 += 4*power(k,-s)*(log(k)-w)**2
print("A1 partial:", nstr(A1, 20))
# tail by fit k*(a ln^2 k + b ln k + c)
xs = []; ys = []
for k in range(K-30, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        w = (a_+b_)/k
        sh += 4*(log(k)-w)**2
    xs.append([k*log(k)**2, k*log(k), k]); ys.append(sh)
A = mpmath.matrix(31, 3); b = mpmath.matrix(31, 1)
for i in range(31):
    for j in range(3): A[i,j] = xs[i][j]
    b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
tail = mpf(0)
for k in range(K+1, 60*K+1):
    tail += power(k,-s)*(c[0,0]*k*log(k)**2 + c[1,0]*k*log(k) + c[2,0]*k)
print("tail:", nstr(tail, 15))
print("A1 full:", nstr(A1+tail, 20))
print("identity value: 0.65682 ; needed for Z'': 0.69018")
'''
exec(code)
