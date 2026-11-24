code = r'''
from squircle import *
from c2 import G_2, G_l, G_x
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
import mpmath
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
s = 4
K = 1200
A1 = A2 = A3 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        lS = log(k)
        A1 += 4*power(k,-s)*(lS-w)**2
        A2 += 4*power(k,-s-2)*L*L
        A3 += 4*power(k,-s-1)*M
# tails via fitted model k*(a ln^2 k + b ln k + c)
def tailfit(vals):
    import mpmath
    A = mpmath.matrix(len(vals), 3); b = mpmath.matrix(len(vals), 1)
    for i, (sh, k) in enumerate(vals):
        A[i,0] = k*log(k)**2; A[i,1] = k*log(k); A[i,2] = k
        b[i,0] = sh
    c = mpmath.lu_solve(A.T*A, A.T*b)
    return c
def shellA1(k):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        sh += 4*(log(k)-w)**2
    return sh*power(k,-s)
def shellA2(k):
    sh = mpf(0)
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_
        sh += 4*L*L
    return sh*power(k,-s-2)
def shellA3(k):
    sh = mpf(0)
    for m0 in range(1, k):
        sh += 4*(m0*log(m0)**2 + (k-m0)*log(k-m0)**2)
    return sh*power(k,-s-1)
for (name, shell, Acur) in [("A1", shellA1, "A1"), ("A2", shellA2, "A2"), ("A3", shellA3, "A3")]:
    vals = [(shell(k), k) for k in range(K-25, K+1)]
    c = tailfit(vals)
    tail = mpf(0)
    for k in range(K+1, 40*K+1):
        tail += c[0,0]*k*log(k)**2*power(k, -s if name=="A1" else (s+2 if name=="A2" else s+1)) + \
                c[1,0]*k*log(k)*power(k, -s if name=="A1" else (s+2 if name=="A2" else s+1)) + \
                c[2,0]*k*power(k, -s if name=="A1" else (s+2 if name=="A2" else s+1))
    print(name, "partial:", nstr(eval(Acur), 20), " tail:", nstr(tail, 12))
    if name=="A1": A1 += tail
    if name=="A2": A2 += tail
    if name=="A3": A3 += tail
print("A1 (s^2/2 T(lnS-w)^2 = s^2/2 * A1?):", nstr(A1, 20))
# careful: A1 partial already has power(k,-s); A1 = sum T (lnS-w)^2 = that
print("U_2 = A2?:", nstr(A2, 20))
print("U_3 = A3?:", nstr(A3, 20))
# A1 partial = sum_k k^{-s} shell(...)/1 : it equals sum T (lnS-w)^2? shell included k^{-s}? yes shellA1 multiplied by power(k,-s)
c2v = (s**2/2)*A1 - (s/2)*(A3 - A2)
print("c_2 = (s^2/2)sumT(lnS-w)^2 - (s/2)(U3-U2) =", nstr(c2v, 20))
print("fit:", 1.573849628813923486848)
'''
exec(code)
