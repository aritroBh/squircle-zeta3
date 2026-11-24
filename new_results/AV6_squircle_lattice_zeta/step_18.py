code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, fsum, power
setdps(60)
p = mpf('0.5'); s = mpf(4); a = s/p
def h(m,n): return power(power(m,p)+power(n,p),-a)
print(type(h(3,7)), type(Jint(mpf(0.4), p, a)), type(Aint(p,a)))
ck = taylor_h(5, 250, p, a, 5)
print([type(c).__name__ for c in ck])
print(type(harmtail(300, s-1)), type(eps_terms(p,s,a,4)[0]))
'''
exec(code)
