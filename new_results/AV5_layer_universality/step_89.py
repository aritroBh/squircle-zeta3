import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr, ellipe
mp.dps = 70
exec(open('step_90.py').read().split('# validate')[0])
Pfull = mpmath.quad(lambda th: h_hex(th, mpf(2)), [0, pi/2, pi, 3*pi/2, 2*pi])
print("P full =", nstr(Pfull,40))
print("ellipse=", nstr(4*ellipe(mpf(2)/3),40))
# check h periodicity
for th in [mpf('0.1'), mpf('0.3')]:
    print(nstr(h_hex(th,mpf(2)),25), nstr(h_hex(th+pi/3,mpf(2)),25), nstr(h_hex(th+pi/6,mpf(2)),25))
