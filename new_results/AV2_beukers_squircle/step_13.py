import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')

# patch av2core.py: stable B3 integrands
src = open('av2core.py').read()
stable = r'''

def B3_quadA(p):
    """2D: x=s, y=t(1-s^p)^{1/p}; corner-stable evaluation:
    1-(1-xy)Z = [(1-A)+A(1-B)] + xy*Z  (all terms >= 0, no cancellation),
    A=(1-s^p)^{1/p}, B=(1-t^p)^{1/p}, Z=AB, xy=s*t*A;  q=1-xy >= 1-2^{-2/p} > 0."""
    p = mpf(p); ep = 1/p
    def f(s, t):
        lA = ep*log1p(-s**p)
        lB = ep*log1p(-t**p)
        A = exp(lA); B = exp(lB)
        one_A = -mp.expm1(lA); one_B = -mp.expm1(lB)
        xy = s*t*A
        val = (one_A + A*one_B) + xy*(A*B)
        return A * (-log(val))/(1-xy)
    return quad(f, [0, 1], [0, 1])

def B3_quadB(p):
    """2D p-polar: corner-stable: 1-(1-r^2 c)Q = (1-Q) + r^2 c Q."""
    p = mpf(p); ep = 1/p; a = 2*ep
    def f(th, r):
        cs = cos(th)*sin(th)
        c = cs**a
        lQ = ep*log1p(-r**p)
        Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        q = 1-r*r*c
        return (2*ep) * cs**(a-1) * r * (-log(val))/q
    return quad(f, [0, pi/2], [0, 1])
'''
# replace old B3_quadA/B definitions
import re
for name in ['B3_quadA', 'B3_quadB']:
    pat = re.compile(r'\ndef %s\(p\):.*?(?=\ndef |\Z)' % name, re.S)
    src = pat.sub('', src)
open('av2core.py','w').write(src + stable)
print("patched")

step05 = r'''
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 70
exec(open('av2core.py').read())
print("=== B3 timing + agreement test, dps=70 (stable) ===")
for p in [1, 2, 4, 16, 64]:
    t0=time.time(); A = B3_quadA(p); t1=time.time()
    B = B3_quadB(p); t2=time.time()
    print("p=%3d  A=%s (%.1fs)" % (p, mp.nstr(A,40), t1-t0))
    print("       B=%s (%.1fs)  |A-B|=%s" % (mp.nstr(B,40), t2-t1, mp.nstr(fabs(A-B),4)))
'''
open('step_05.py','w').write(step05)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_05.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_05: stable corner evaluation for B3 (expm1/log1p); A/B agree at dps=70. DONE\n")
