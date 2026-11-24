code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
mpows = [power(m, p) for m in range(1, 3000)]
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-80); m = 0
    while True:
        m += 1
        e = exp(-t*mpows[m-1])
        f += e
        if e < cutoff: break
    return 1+2*f
def F(t):
    th = theta(t); return th*th-1
delta = mpf('2.46e-150'); ld = log(delta)
# quad with different maxdegree / splitting
I1a = quad(lambda x: exp(a*x)*F(exp(x)), [ld, 0])
I1b = quad(lambda x: exp(a*x)*F(exp(x)), [ld, -100, -20, 0], maxdegree=12)
I1c = quad(lambda x: exp(a*x)*F(exp(x)), [ld, -50, 0], maxdegree=10)
print("I1a:", nstr(I1a, 40))
print("I1b:", nstr(I1b, 40))
print("I1c:", nstr(I1c, 40))
I2 = quad(lambda t: t**(a-1)*F(t), [1, inf])
print("I2 :", nstr(I2, 40))
'''
exec(code)
