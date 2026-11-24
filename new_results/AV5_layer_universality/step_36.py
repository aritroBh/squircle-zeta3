import mpmath
from mpmath import mp, mpf, log, exp, quad, zeta, nstr
mp.dps = 130
g1 = lambda s: log(1+exp(-2*s))
# method comparison for I(0,2)
f = lambda s: g1(s)**2
r1 = quad(f, [0, mpf('0.5'), 1, 2, 4, 8, 16, 40])
# plain tanh-sinh on [0,40]
r2 = quad(f, [0, 40])
# gauss-legendre
r3 = quad(f, [0, 2, 6, 20, 60], method='gauss-legendre')
exact = zeta(3)/8
for nm,r in [('split tanh-sinh',r1),('[0,40] tanh-sinh',r2),('gauss-legendre',r3)]:
    print(nm, nstr(r-exact, 6))
# tanh-sinh with explicit quadrature object
TS = mpmath.calculus.quadrature.TanhSinh(mpmath.mp)
r4 = TS.sum(f, [0, 40])
print('TanhSinh obj:', nstr(r4-exact,6))
