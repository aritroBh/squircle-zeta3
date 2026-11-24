import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb, acb
ctx.dps = 40
from step_07_partb import *
for p in [0, 1, 2, 3, 4, 6, 9]:
    v = mu_p_of(400, p)
    print(f"mu({p}) = {v.str(8)}")
ms = momentsum_of(400)
print("momentsum =", ms.str(8))
print("ModelRatio =", ModelRatio_of(400).str(10))
print("T1ratio =", T1ratio_of(400).str(10))
'''
open("t_mu.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_mu.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")

from mpmath import mp, mpf, pi, sqrt, exp, beta
mp.dps = 50
j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
def df2(k):
    r = mpf(1)
    while k > 1: r *= k; k -= 2
    return r
N = mpf(400); y = sqrt(pi/(2*N))
def ell(m):
    if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
    return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
l1,l2 = (ell(m) for m in range(1,3))
t0 = -l1/(2*l2)
def mu(p):
    base = exp(-l2*t0**2) * sqrt(-pi/l2)
    return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
for p in [0,1,2,3,4,6,9]:
    print(f"mpmath mu({p}) = {mp.nstr(mu(p), 10)}")
