code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma, fsum, hyper
setdps(70)
p = mpf(50); a = mpf(4)/p
for x in [mpf(41)/250*250/41, mpf('6.09756'), mpf('48.7805'), mpf('0.5'), mpf('1.0'), mpf('3.0')]:
    j1 = Jint(x, p, a)
    j2 = quad(lambda u: (1+u**p)**(-a), [x, inf])
    print(f"x={nstr(x,6)}: Jint={nstr(j1,30)}")
    print(f"          quad={nstr(j2,30)}  diff={nstr(j1-j2,4)}")
print("Aint:", nstr(Aint(p,a), 30))
print("quad:", nstr(quad(lambda u: (1+u**p)**(-a), [0, inf]), 30))
'''
exec(code)
