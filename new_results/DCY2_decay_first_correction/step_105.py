import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_sweep.py").read()
old = '''    def GsplitB(jj, gam_lo):
        # On [0,T], T = c2m/(2 gam): e^{gam|t|^3} <= e^{(c2m/2) t^2} -> widened moment.
        # On [T,Delta]: the damped-Gaussian/polynomial edge.
        gamD_lo = gam_lo*D_lo
        T_lo = min(D_lo, c_lo/(2*gam_lo))
        base = Mtilde_of(Nhi, jj, c_hi/2)
        if gamD_lo <= c_lo/2:
            return base
        excess_hi = gamD_lo - c_hi
        e2 = (ReL1max_of(Nlo)*D_lo).exp()
        if excess_hi < 0:
            rate_lo = c_hi - gamD_lo
            edge = D_hi**(jj+1)*(-rate_lo*T_lo**2).exp()*e2/(jj+1)
        else:
            edge = e2*(excess_hi*D_hi**2).exp()*D_hi**(jj+1)/(jj+1)
        return base + edge'''
new = '''    def GsplitB(jj, gam_lo):
        # base: min( e^{gam T^3} Mtilde(c2m), widened Mtilde(c2m/2) ) on [0,T];
        # edge [T,Delta]: damped incomplete-Gamma integral (excess<0) or
        # polynomial moment with the excess factor (excess>0).
        gamD_lo = gam_lo*D_lo
        T_lo = min(D_lo, c_lo/(2*gam_lo))
        base = min((gam_lo*T_lo**3).exp()*Mtilde_of(Nhi, jj),
                   Mtilde_of(Nhi, jj, c_hi/2))
        if gamD_lo <= c_lo/2:
            return base
        excess_hi = gamD_lo - c_hi
        e2 = (ReL1max_of(Nlo)*D_lo).exp()
        if excess_hi < 0:
            rate_lo = c_hi - gamD_lo
            edge = EdgeInt_of(Nlo, jj, rate_lo, T_lo)*e2
        else:
            edge = e2*(excess_hi*D_hi**2).exp()*D_hi**(jj+1)/(jj+1)
        return base + edge'''
assert old in src
src = src.replace(old, new)
open("dcy2_sweep.py","w").write(src)
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
open("t_sweep3.py","w").write(sweeprun)
r = subprocess.run(["<project>/.venv/bin/python","t_sweep3.py"],
                   capture_output=True, text=True, cwd=".", timeout=1200)
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
