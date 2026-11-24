code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma, gammainc, fsum
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
ap = gamma(1+1/p); ga = gamma(a)
delta = mpf('2.46e-150')
cache = {}
def mpw(m):
    v = cache.get(m)
    if v is None:
        v = power(m, p); cache[m] = v
    return v
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpw(m)); f += e
        if e < cutoff: break
    return 1+2*f
def F(t):
    th = theta(t); return th*th-1
q1 = quad(lambda u: u**(a-1)*F(delta*u), [0, 1]) * delta**a
small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
print("small quad:", nstr(q1, 40))
print("small form:", nstr(small, 40))
print("diff:", nstr(q1-small, 6))
'''
exec(code)
