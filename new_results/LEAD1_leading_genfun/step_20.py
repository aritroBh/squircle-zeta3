import time, mpmath
mpmath.mp.dps = 50
from mpmath import mpf, sqrt, log, quad, pi
a = 1 + sqrt(2); b = sqrt(2) - 1
def tfun(v): return (v*v + 2*v - 1)/(4*v)
def integrand(v, i, j, r):
    t = tfun(v)
    return log(t)**i * log(1-t)**j / (v - r)
t0=time.time()
N = {}
for r,nm in [(-a,'-a'),(mpf(0),'0'),(-b,'-b'),(a,'a')]:
    N[nm] = quad(lambda v: integrand(v,1,0,r), [b,1,a])
c2 = -N['-a'] + N['0']/2 - N['-b']/2 + N['a']/2
lam = log(2)
print('c2      =', c2)
print('target  =', (pi**2 + 6*lam**2)/24)
print('diff    =', c2 - (pi**2 + 6*lam**2)/24)
print(f'total {time.time()-t0:.2f}s')
