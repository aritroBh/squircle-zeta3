import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
from mpmath import mp, mpf, log, polylog, fabs, nstr, loggamma, exp, quad
mp.dps = 80
exec(open('av2core.py').read())
p = mpf(1024)
FP = FPrecomp(p, '0.02', 500)
def g(c):
    return (polylog(2, c) + log(c)*mp.log1p(-c))/c
def H(c):
    return FP.F(c) - g(c)/2
cstar = 2**(-2/p)
E_direct = 2*quad(lambda c: H(c)/mp.sqrt(1-4*c**p) if 1-4*c**p > 0 else 0, [0, cstar])
print("E from integral:", nstr(E_direct, 20))
B3 = mpf('2.40408689560606935868850025202647004088848179266562375168724')
B2 = B2_quad2(1024)
print("E measured     :", nstr(B3-2*B2, 20))
# examine H at small c
for c in [1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1]:
    c = mpf(c)
    print("c=%8.0e  H=%s  (1/2)log c + p/4 = %s" % (c, nstr(H(c),12), nstr(log(c)/2+p/4, 12)))
