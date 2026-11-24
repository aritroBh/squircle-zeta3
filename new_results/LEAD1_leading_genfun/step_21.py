import time, mpmath
mpmath.mp.dps = 400
from mpmath import mpf, sqrt, log, quad, pi
a = 1 + sqrt(2); b = sqrt(2) - 1
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v)
    u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
t0=time.time()
v1 = quad(lambda v: integrand(v, 8, 0, -a), [b, 1, a])
t1=time.time()
v2 = quad(lambda v: integrand(v, 4, 4, mpf(0)), [b, 1, a])
t2=time.time()
print('N_(8,0),-a =', mpmath.nstr(v1, 30), f'{t1-t0:.2f}s')
print('N_(4,4),0  =', mpmath.nstr(v2, 30), f'{t2-t1:.2f}s')
