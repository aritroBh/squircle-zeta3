import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb, acb
ctx.dps = 40
from step_07_partb import *
print("S_of(400) =", S_of(400).str(10))
print("exp(S) =", S_of(400).exp().str(10))
print("H_of(400) =", H_of(400).str(10))
print("momentsum =", momentsum_of(400).str(10))
'''
open("t_factors.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_factors.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")

from mpmath import mp, mpf, pi, sqrt, exp, beta
mp.dps = 50
j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
def Hc(s): return beta(s/2, mpf(3)/2)/2
N = mpf(400); y = sqrt(pi/(2*N)); s0 = omega/y
S = sum((-1)**m*pi*y**(m-2)*wb**m/(2*m) for m in range(3,12))
print("mpmath S =", mp.nstr(S, 12))
print("mpmath exp(S) =", mp.nstr(exp(S), 12))
print("mpmath H0 =", mp.nstr(Hc(s0), 12))
