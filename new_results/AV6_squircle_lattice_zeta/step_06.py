code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, diff, fsum, power
setdps(50)
p = mpf(2); s = mpf(4); a = s/p

# check taylor_h against numerical derivative
m, y = 5, 120
ck = taylor_h(m, y, p, a, 6)
g = lambda t: (m**p + (y+t)**p)**(-a)
for k in range(5):
    exact = diff(g, 0, n=k)/fac(k)
    print(k, nstr(ck[k], 12), nstr(exact, 12), nstr(fabs(ck[k]-exact)/fabs(exact), 3))
'''
exec(code)
