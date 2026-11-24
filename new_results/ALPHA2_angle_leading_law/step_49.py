import sys, os
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from mpmath import mp, mpf, nstr, pi, sin, cos, acos, gamma, hyper, diffs, factorial, polylog
mp.dps = 60
# assembled values saved from previous run? recompute quickly for c12 and store
from assemble import assemble_cw
cw, c, r = assemble_cw('c12', wmax=7)
# T(alpha) check
c2 = (mpf(c.numerator)/c.denominator)**2
T = pi**2/12 - polylog(2, c2)/2
print("c(2) assembled =", nstr(cw[2], 40))
print("T(alpha)       =", nstr(T, 40))
print("|diff| =", nstr(abs(cw[2]-T), 5))
