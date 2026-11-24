import sys, os
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1); os.chdir(fd1)
from mpmath import mp, mpf, nstr, pi, acos
mp.dps=50
from assemble import assemble_cw
import corner_series
cw, c, r = assemble_cw('c12', wmax=7, dps=50)
print("c,r:", c, r)
a = 2*acos(mpf(c.numerator)/c.denominator)
M = corner_series.corner_coeffs(a, 6)
for w in range(2,8):
    print(f"w={w}: c(w)={nstr(cw[w],25)}  M_{w-1}={nstr(M[w-1],25)}  ratio={nstr(cw[w]/M[w-1],8)}")
