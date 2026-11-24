import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step34 = r'''
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, euler, nsum, inf
mp.dps = 105
z2, z3, ln2, ln3 = zeta(2), zeta(3), log(2), log(3)
from math import comb
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag
def Hm(m):
    return nsum(lambda k: 1/k, [1, int(m)])
S1  = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * Hm(m)/(2*m+1), [1, inf])
S1p = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * Hm(2*m)/(2*m+1), [1, inf])
S2  = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m)/(2*m+1)**2, [1, inf])
S3  = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m)/(2*m+1), [1, inf])
print("S1  =", nstr(S1, 55))
print("S1' =", nstr(S1p, 55))
print("S2  =", nstr(S2, 30), " vs Cl2-1:", nstr(S2-(Cl2-1), 10))
print("S3  =", nstr(S3, 30), " vs (pi-3)/3:", nstr(S3-(pi-3)/3, 10))
def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**10):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    ok = rel is not None and rel[0] != 0
    print("%s: %s" % (name, rel if ok else "(none)"))
hunt("S1", S1, [mpf(1), ln2, ln2**2, ln3, pi, pi*ln2, Cl2, z3])
hunt("S1b", S1, [mpf(1), pi*Cl2, Cl2, sqrt(3)*Cl2, ln3, ln2, pi])
hunt("S1p", S1p, [mpf(1), ln2, ln2**2, ln3, pi, pi*ln2, Cl2, z3])
hunt("S1pb", S1p, [mpf(1), pi*Cl2, Cl2, sqrt(3)*Cl2, ln3, ln2, pi])
exec(open('av2core.py').read())
B2_2, _ = B2_series(2)
E2 = B2_2 - (mpf(3)/32 + z3/12)
print("E2 =", nstr(E2, 55))
E2check = pi/8*(S1 - S1p + S2 + 2*ln2*S3)
print("E2 vs (pi/8)(S1-S1'+S2+2ln2 S3):", nstr(E2 - E2check, 10))
hunt("E2a", E2, [mpf(1), pi, pi*ln2, pi*ln2**2, pi**2, Cl2, pi*Cl2, ln3])
hunt("E2b", E2/pi, [mpf(1), ln2, ln2**2, pi, Cl2, ln3])
'''
open('step_34.py','w').write(step34)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_34.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_34: S2=Cl2(pi/3)-1, S3=(pi-3)/3; hunts for S1, E2. DONE\n")
