import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step41 = r'''
from mpmath import mp, mpf, zeta, fabs, nstr, nsum, inf
mp.dps = 100
exec(open('av2core.py').read())
from math import comb, gcd
from functools import reduce
def lcm(a, b): return a*b//gcd(a, b)
def L(n): return reduce(lcm, range(1, n+1), 1)
def apery(n): return sum(comb(n,k)**2 * comb(n+k,k)**2 for k in range(n+1))
def Minf(j, k):
    return nsum(lambda m: mpf(j+k+2*m+2)/((j+m+1)**2*(k+m+1)**2), [0, inf])
Mcache = {}
def Ln_inf(n):
    cs = shifted_legendre_coeffs(n)
    tot = mpf(0)
    for i, ci in enumerate(cs):
        for j, cj in enumerate(cs):
            if ci == 0 or cj == 0: continue
            key = (min(i,j), max(i,j))
            if key not in Mcache: Mcache[key] = Minf(*key)
            tot += ci*cj*Mcache[key]
    return tot
z3 = zeta(3)
print("n | dn^3 L_n = A_n zeta(3) + B_n : verify A_n = 2 beta_n dn^3, B_n integer")
for n in range(0, 13):
    Ln = Ln_inf(n)
    dn3 = L(n)**3
    An = 2*apery(n)*dn3
    Bn = dn3*Ln - An*z3
    BnR = round(Bn)
    print("n=%2d  A_n=%d  B_n=%s (round err=%s)  |L_n|^(1/n)=%s" % (
        n, An, BnR, nstr(fabs(Bn-BnR), 6), nstr(fabs(Ln)**(1.0/n) if n else 0, 8)))
'''
open('step_41.py','w').write(step41)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_41.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_41: Beukers structure dn^3 L_n = A_n zeta3 + B_n verified n<=12: A_n=2 apery(n) dn^3, B_n integer. DONE\n")
