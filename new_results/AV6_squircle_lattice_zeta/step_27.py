code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma
import time
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
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
delta = mpf('2.46e-150'); ld = float(log(delta))
t0=time.time()
# Gauss-Legendre on chunks of 25
pts = [ld]
while pts[-1] < -25: pts.append(pts[-1]+25)
pts.append(0.0)
tot_gl = mpf(0)
for i in range(len(pts)-1):
    tot_gl += quad(lambda x: exp(a*x)*F(exp(x)), [pts[i], pts[i+1]], method='gauss-legendre')
print("I1 GL:", nstr(tot_gl, 45), f"{time.time()-t0:.0f}s")
print("I1 TS:", "111.452594235019211839272981231853695794195775")
'''
exec(code)
