import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")

# R4a: E_tailmom x2 (two-sided)
src = open("dcy2_sweep.py").read()
old = """    pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
          6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
    Et = arb(0)
    for p, cf in pp.items():
        Et += cf*TailMom_of(Nhi, p)
    return E + Et"""
new = """    pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
          6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
    Et = arb(0)
    for p, cf in pp.items():
        Et += cf*TailMom_of(Nhi, p)
    return E + 2*Et   # two-sided window tail (R4)"""
assert old in src
src = src.replace(old, new)
open("dcy2_sweep.py","w").write(src)

# R4b: Mtilde correction 2 e^{a^2/4c} t_R^{j+1}/(j+1)
src = open("dcy2_cell.py").read()
old = """def Mtilde_of(N, jj, c2use=None):
    c2m = c2min_of(N) if c2use is None else c2use
    ReL1 = ReL1max_of(N); tR = ReL1/(2*c2m)
    tot = arb(0)
    for i in range(jj+1):
        tot += binom(jj,i)*tR**(jj-i)*c2m**(-(arb(i)+1)/2)*((arb(i)+1)/2).gamma()
    return (ReL1**2/(4*c2m)).exp()*tot"""
new = """def Mtilde_of(N, jj, c2use=None):
    c2m = c2min_of(N) if c2use is None else c2use
    ReL1 = ReL1max_of(N); tR = ReL1/(2*c2m)
    tot = arb(0)
    for i in range(jj+1):
        tot += binom(jj,i)*tR**(jj-i)*c2m**(-(arb(i)+1)/2)*((arb(i)+1)/2).gamma()
    # R4: provable correction for the [-t_R, 0] edge of the shifted Gaussian
    corr = 2*(ReL1**2/(4*c2m)).exp()*tR**(jj+1)/(jj+1)
    return (ReL1**2/(4*c2m)).exp()*tot + corr"""
assert old in src
src = src.replace(old, new)
open("dcy2_cell.py","w").write(src)
print("R4 applied")

sweeprun = '''
import time
from flint import ctx, arb
ctx.dps = 40
from dcy2_sweep import LHS_cell
t0 = time.monotonic()
worst = arb(0); worstcell = None
x = 400.0
while x < 100000.0:
    x2 = min(100000.0, x*1.02)
    v = LHS_cell(x, x2)
    if v.upper() > worst:
        worst = v.upper(); worstcell = (x, x2)
    x = x2
print(f"worst LHS upper = {worst.str(10)} at {worstcell}")
print(f"margin vs C2=1e6: {float(10**6/worst):.3f}x")
print(f"elapsed: {time.monotonic()-t0:.1f}s")
'''
open("t_sweep_r4.py","w").write(sweeprun)
r = subprocess.run(["<project>/.venv/bin/python","t_sweep_r4.py"],
                   capture_output=True, text=True, cwd=".", timeout=1200)
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
