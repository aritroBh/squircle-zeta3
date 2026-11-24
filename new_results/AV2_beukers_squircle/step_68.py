import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step32 = r'''
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, euler, diff
mp.dps = 110
exec(open('av2core.py').read())
B2_1, _ = B2_series(1)
B3_1 = mpf('0.24661730238492374056516806883959205599718610151827')
B2_2, _ = B2_series(2)
z2, z3, z4, z5 = zeta(2), zeta(3), zeta(4), zeta(5)
ln2, ln3, gam = log(2), log(3), euler
zp2 = diff(zeta, 2)   # zeta'(2)
zp3 = diff(zeta, 3)
Li3h = polylog(3, mpf(1)/2)
print("-zeta'(2) =", nstr(-zp2, 30), "  -zeta'(3) =", nstr(-zp3, 25))
def hunt(name, val, basis, tol=mpf('1e-55'), maxcoeff=10**11):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    ok = rel is not None and rel[0] != 0
    print("%s: %s" % (name, rel if ok else "(none)"))
print("=== B3(1) with zeta' ===")
hunt("b31", B3_1, [mpf(1), z2, z3, zp2, zp3])
hunt("b31b", B3_1, [mpf(1), z2, z3, zp2, zp3, ln2, z2*ln2, z3*ln2, ln2**2, ln2**3])
hunt("b31c", B3_1, [mpf(1), z2, z3, zp2, ln2, z2*ln2, ln2**2, Li3h])
hunt("b31d", B3_1, [mpf(1), z3, z2, gam, ln2, ln3])
print("=== B2(2) with ln2 family ===")
hunt("b22a", B2_2, [mpf(1), ln2, ln2**2, ln2**3, pi, pi*ln2, pi*ln2**2, z2, z3, z2*ln2])
hunt("b22b", B2_2, [mpf(1), ln2, ln2**2, pi, pi*ln2, z2, z3, z2*ln2, Li3h, zp2])
hunt("b22c", B2_2, [mpf(1), ln2, pi, pi*ln2, z3, zp2, ln2**2, ln2**3])
print("=== sanity: B2(1) = 2 zeta(3)/3 at 110 digits ===")
print(nstr(B2_1 - 2*z3/3, 30))
'''
open('step_32.py','w').write(step32)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_32.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_32: richer PSLQ hunts (zeta' etc). DONE\n")
