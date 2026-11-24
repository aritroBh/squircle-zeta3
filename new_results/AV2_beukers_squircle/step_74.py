import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step37 = r'''
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, nsum, inf
mp.dps = 110
from math import comb
z2, z3, ln2, ln3 = zeta(2), zeta(3), log(2), log(3)
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag
def H(m):
    return nsum(lambda k: 1/k, [1, int(m)])
S1  = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * H(m)/(2*m+1), [1, inf])
S1p = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * H(2*m)/(2*m+1), [1, inf])
def hunt(name, val, basis, tol=mpf('1e-50'), maxcoeff=10**8):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if rel is None or rel[0] == 0:
        print("%s: (none)" % name); return
    resid = fabs(sum(mpf(rel[i])*([val]+basis)[i] for i in range(len(rel))))
    scale = max(fabs(mpf(rel[i])*([val]+basis)[i]) for i in range(len(rel)))
    tag = "OK" if resid/scale < mpf('1e-40') else "WEAK"
    print("%s [%s]: %s  resid=%s" % (name, tag, rel, nstr(resid/scale, 5)))
big = [mpf(1), pi, pi**2, ln2, ln3, ln2**2, Cl2, pi*Cl2, sqrt(3)*Cl2, pi*ln3, z3]
hunt("S1", S1, big)
hunt("S1'", S1p, big)
small = [mpf(1), pi, ln2, Cl2, pi*Cl2, z3]
hunt("S1 s", S1, small)
hunt("S1' s", S1p, small)
# guess structure: S1 = a + b pi + c Cl2 + d pi Cl2?  print key combos
print("S1 =", nstr(S1, 50))
print("S1'=", nstr(S1p, 50))
print("2 S1 - S1' =", nstr(2*S1-S1p, 40))
print("S1+S1' =", nstr(S1+S1p, 40))
print("S1'-S1 =", nstr(S1p-S1, 40))
'''
open('step_37.py','w').write(step37)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_37.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
