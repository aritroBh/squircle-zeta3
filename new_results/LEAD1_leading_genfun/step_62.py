import mpmath
from mpmath import mpf, sqrt, log, quad
a = 1+sqrt(2); b = sqrt(2)-1
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v); u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
for dps in [100, 200, 300, 400]:
    mpmath.mp.dps = dps
    aa = 1+sqrt(2); bb = sqrt(2)-1
    v1 = quad(lambda v: integrand(v, 9, 0, -aa), [bb,1,aa])
    v2 = quad(lambda v: integrand(v, 5, 4, mpf(0)), [bb,1,aa])
    print(dps, mpmath.nstr(v1, 25), mpmath.nstr(v2, 25))
