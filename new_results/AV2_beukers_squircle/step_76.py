import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step39 = r'''
from mpmath import mp, mpf, zeta, fabs, log, nstr, pslq, loggamma, digamma, exp, euler
mp.dps = 110
z2, z3, ln2 = zeta(2), zeta(3), log(2)

def dB2_deps(eps):
    """d/d eps of B2 = sum T_n, T_n = eps^3 G(n eps)^2/G(2n eps+1) (psi(2n eps+1)-psi(n eps))."""
    S = mpf(0); n = 1
    tol = mpf(10)**(-(mp.dps-10))
    while True:
        w = n*eps
        lg = 2*loggamma(w) - loggamma(2*w+1)
        T = eps**3*exp(lg)*(digamma(2*w+1)-digamma(w))
        # dT/deps = T*(3/eps + 2n psi(w) - 2n psi(2w+1)) + eps^3 e^lg (2n psi'(2w+1)- n psi'(w))
        from mpmath import psi
        dT = T*(3/eps + 2*n*digamma(w) - 2*n*digamma(2*w+1)) + eps**3*exp(lg)*(2*n*psi(1, 2*w+1) - n*psi(1, w))
        S += dT
        if n > 40/eps and fabs(dT) < tol:
            break
        n += 1
    return S

d = dB2_deps(mpf(1))
print("dB2/deps at eps=1 =", nstr(d, 60))
print("d1 = -dB2/deps =", nstr(-d, 60))
def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**9):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if rel is None or rel[0] == 0:
        print("%s: (none)" % name); return
    resid = fabs(sum(mpf(rel[i])*([val]+basis)[i] for i in range(len(rel))))
    scale = max(fabs(mpf(rel[i])*([val]+basis)[i]) for i in range(len(rel)))
    tag = "OK" if resid/scale < mpf('1e-40') else "WEAK"
    print("%s [%s]: %s  resid=%s" % (name, tag, rel, nstr(resid/scale, 5)))
hunt("d1", -d, [mpf(1), z2, z3, z2*z3, ln2, z2*ln2, z3*ln2, ln2**2, euler, euler*z2])
hunt("d1b", -d, [mpf(1), z2, z3, zeta(4), z2*log(2), zeta(4)*log(2), log(3), z2*log(3)])
# also fit d2 numerically from the tau series via exact second derivative? skip; hunt d2 from tau fit
'''
open('step_39.py','w').write(step39)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_39.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_39: exact dB2/deps at eps=1; d1 hunts. DONE\n")
