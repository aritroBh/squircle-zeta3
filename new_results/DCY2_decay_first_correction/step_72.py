import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import arb, ctx
ctx.dps = 40
PI = arb.pi()
def Gamma_upper_half(n, x):
    res = PI.sqrt()*(x.sqrt()).erfc()
    G = res
    for k in range(1, n+1):
        G = (k-arb(1)/2)*G + x**(k-arb(1)/2)*(-x).exp()
    return G
def Gamma_upper_int(n, x):
    tot = arb(0)
    for i in range(n):
        tot += x**i/arb(i).fac()
    return arb(n-1).fac()*(-x).exp()*tot
for s2, x in [(13, "0.742"), (13, "3.08"), (3, "1.0"), (9, "4.757"), (1, "0.5"), (13, "0.28"), (7, "12.0")]:
    n = (s2-1)//2
    a = arb(x)
    Gh = Gamma_upper_half(n, a)
    Gi = Gamma_upper_int(n, a) if n+1 == (s2+1)//2 and (s2 % 2 == 1) else None
    print(f"Gamma({s2}/2, {x}): half-int impl = {Gh.str(12)}")
for n, x in [(7, "0.742"), (4, "3.08"), (2, "1.0"), (1, "0.5")]:
    a = arb(x)
    print(f"Gamma({n}, {x}): int impl = {Gamma_upper_int(n, a).str(12)}")
'''
open("t_gamma.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_gamma.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")

from mpmath import mp, gammainc, mpf
mp.dps = 40
for s, x in [(6.5, 0.742), (6.5, 3.08), (1.5, 1.0), (4.5, 4.757), (0.5, 0.5), (6.5, 0.28), (3.5, 12.0)]:
    print(f"mpmath Gamma({s}, {x}) =", gammainc(s, x, regularized=False))
for n, x in [(7, 0.742), (4, 3.08), (2, 1.0), (1, 0.5)]:
    print(f"mpmath Gamma({n}, {x}) =", gammainc(n, x, regularized=False))
