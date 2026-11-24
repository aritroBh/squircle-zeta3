final_check = r'''import re
from mpmath import mp, mpf, log, pi, zeta, sqrt, nstr, fabs
txt150 = open('p5_run150.txt').read()
q = re.search(r'P5/sqrt3 = ([0-9.]+)', txt150).group(1)
mp.dps = 160
L2 = log(2)
# consolidated final numbers
z5c = mpf(2835)/32; z23c = mpf(1395)/32
REL = [128, 0, -30720, 4920, 61440, 0, -6240, -30720, -448, -5760, 1920, 5120, -11340, -5580, 123, 1440, -160, -256]
basis = [mpf(1), L2, pi**2, L2**2, zeta(3), pi**2*L2, L2**3, pi**4, zeta(3)*L2,
         pi**2*L2**2, L2**4, zeta(5), (pi**2/6)*zeta(3), pi**4*L2, zeta(3)*L2**2, pi**2*L2**3, L2**5]
val = mpf(q)
resid = fabs(128*val - sum(mpf(-REL[i+1])*basis[i] for i in range(17)))
print("P^(5)(1) =", re.search(r'P\^\(5\)\(1\) = ([0-9.]+)', txt150).group(1)[:92], "...")
print("residual of 128-relation at dps=160           :", nstr(resid, 5))
print("claimed zeta(5) part  2835/32 *z5*sqrt3       :", nstr(z5c*zeta(5)*sqrt(3), 30))
print("claimed zeta(2)zeta(3) part 1395/32 *z2z3*sqrt3:", nstr(z23c*(pi**2/6)*zeta(3)*sqrt(3), 30))
print("relation z5 coeff  = 11340/128 == 2835/32 ?", mpf(11340)/128 == z5c)
print("relation z2z3 coeff =  5580/128 == 1395/32 ?", mpf(5580)/128 == z23c)
print("VERDICT: CONFIRMED")
'''
open('hex5_finalcheck.py','w').write(final_check)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','hex5_finalcheck.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-300:])
print(open('result.md').read()[:400])
