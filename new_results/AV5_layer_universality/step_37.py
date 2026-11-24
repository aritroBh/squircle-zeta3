import mpmath, pickle
from mpmath import mp, mpf, log, exp, quad, zeta, nstr
mp.dps = 130
g1 = lambda s: log(1+exp(-2*s))
f = lambda s: g1(s)**2
exact = zeta(3)/8
for nsplit in [10, 20, 40, 80]:
    pts = [mpf(40)*j/nsplit for j in range(nsplit+1)]
    r = quad(f, pts, method='gauss-legendre')
    print(nsplit, nstr(r-exact,6))
