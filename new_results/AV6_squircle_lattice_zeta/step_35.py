code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma
import time
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
ap = gamma(1+1/p); ga = gamma(a)
print("ga = Gamma(0.08) =", nstr(ga, 30))
ZA_ = ZA(50, 4)
print("ZA * ga =", nstr(ZA_*ga, 45))
delta = mpf('2.46e-150'); ld = float(log(delta))
mpows = [power(m, p) for m in range(1, 5000)]
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpows[m-1]); f += e
        if e < cutoff: break
    return 1+2*f
def F(t):
    th = theta(t); return th*th-1
small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
print("small =", nstr(small, 30))
t0=time.time()
pts = [ld]
while pts[-1] < -25: pts.append(pts[-1]+25)
pts.append(0.0)
I1 = mpf(0)
for i in range(len(pts)-1):
    I1 += quad(lambda x: exp(a*x)*F(exp(x)), [pts[i], pts[i+1]], method='gauss-legendre')
print("I1 (GL chunks) =", nstr(I1, 45), time.time()-t0)
I2x = quad(lambda x: exp(a*x)*F(exp(x)), [0.0, 5.0])   # = int_1^inf t^{a-1} F(t) dt
print("I2 (x-form) =", nstr(I2x, 45))
total = small + I1 + I2x
print("total =", nstr(total, 45))
print("total/ga =", nstr(total/ga, 45))
print("ZA       =", nstr(ZA_, 45))
'''
exec(code)
