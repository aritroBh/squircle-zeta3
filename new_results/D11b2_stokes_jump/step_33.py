# ======== STEP 2c: growth/asymptotics of f_n ========
from mpmath import mp, mpf, pi, factorial
import json
mp.dps = 120
F = [mpf(v) for v in json.load(open('step2_num200.json'))['f200']]
def pat(n):   # predicted sign pattern sqrt2 cos(pi(3-2n)/4) = +-1: n%4: 0:-, 1:+, 2:+, 3:-
    return [ -1, 1, 1, -1 ][n % 4]
print(" n    R_n = f_n/(n!(2/pi)^n)          C_n = R_n*n^{3/2}     C_n/pattern  (even n)")
for n in range(10, 61, 2):
    R = F[n]/(factorial(n)*(2/pi)**n)
    C = R*n**mpf('1.5')
    print(f"{n:3d}  {mp.nstr(R, 12):>18}  {mp.nstr(C, 12):>14}   {mp.nstr(C/pat(n), 12)}")
print("\n (odd n)")
for n in range(11, 60, 2):
    R = F[n]/(factorial(n)*(2/pi)**n)
    C = R*n**mpf('1.5')
    print(f"{n:3d}  {mp.nstr(R, 12):>18}  {mp.nstr(C, 12):>14}   {mp.nstr(C/pat(n), 12)}")
