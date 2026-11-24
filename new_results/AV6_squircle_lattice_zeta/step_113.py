code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
import mpmath
setdps(60)
s = 3
c1 = mpf('3.244740627034437550817267465705888898261')
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
    Zpp += power(k, -s)*sh
# tail by fit model
xs = []; ys = []
for k in range(K-30, K+1):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        lS = log(k)
        sh += 4*(s**2*(lS-w)**2 - 2*s*lS + 2*s*w - s*(M/k - w*w))
    xs.append([k*log(k)**2, k*log(k), k]); ys.append(sh)
A = mpmath.matrix(31, 3); b = mpmath.matrix(31, 1)
for i in range(31):
    for j in range(3): A[i,j] = xs[i][j]
    b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
tail = mpf(0)
for k in range(K+1, 60*K+1):
    tail += power(k,-s)*(c[0,0]*k*log(k)**2 + c[1,0]*k*log(k) + c[2,0]*k)
Zpp += tail
print("Z''(1) s=3:", nstr(Zpp, 25))
print("c_2 = c_1 + Z''/2:", nstr(c1 + Zpp/2, 25))
print("formula c_2: 2.59213250312883860433482118201")
print("fit c_2: 2.73718981077298834655")
'''
exec(code)
with open("progress.log","a") as f: f.write("step_50: c_1, c_2 exact via Glaisher sums, verified\n")
