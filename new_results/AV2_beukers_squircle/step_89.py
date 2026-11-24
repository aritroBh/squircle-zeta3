import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step46 = r'''
from mpmath import mp, mpf, zeta, fabs, nstr
mp.dps = 60
from fractions import Fraction
from math import comb, factorial

# dB3/dp at p=1 via moments:
# J(j,j,m,1) = j!^2 m!/(2j+m+3)!
# dJ/dp|_1 = J * [-3 - 2(j+1)H_j - (m+1)H_m + (2j+m+3)H_{2j+m+3}]
# dB3/dp = sum_m sum_j (-1)^j C(m,j) dJ/dp
def H(n):
    return Fraction(1, 1) if n == 1 else sum(Fraction(1, k) for k in range(1, n+1))
N = 220
Hs = [Fraction(0)] + [H(n) for n in range(1, 2*N+60)]
def Jval(j, m):
    return Fraction(factorial(j)**2 * factorial(m), factorial(2*j+m+3))
S = Fraction(0)
Km = []
for m in range(0, N):
    sm = Fraction(0)
    for j in range(0, m+1):
        J = Jval(j, m)
        dJ = J * (Fraction(-3) - 2*(j+1)*Hs[j] - (m+1)*Hs[m] + (2*j+m+3)*Hs[2*j+m+3])
        sm += Fraction((-1)**j * comb(m, j)) * dJ
    S += sm
    Km.append(sm)
    if m % 40 == 0:
        print("m=%d partial=%s K_m=%s" % (m, nstr(mpf(S.numerator)/mpf(S.denominator), 20),
              nstr(mpf(sm.numerator)/mpf(sm.denominator), 8)))
print("partial sum (m<%d) =" % N, nstr(mpf(S.numerator)/mpf(S.denominator), 40))
# tail: K_m ~ (ln m)^2/m^2-ish; fit last chunk: log K_m vs log m
import math
for m in [N-40, N-30, N-20, N-10, N-5, N-1]:
    v = Km[m]
    print("K_%d = %s  m^2 K_m = %s" % (m, nstr(mpf(v.numerator)/mpf(v.denominator), 10),
          nstr(mpf(v.numerator)/mpf(v.denominator)*m*m, 10)))
'''
open('step_46.py','w').write(step46)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_46.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
