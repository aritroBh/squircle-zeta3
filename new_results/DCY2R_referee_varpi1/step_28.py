import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 45
t0 = time.monotonic()
from dcy2_sweep import LHS_cell
C2 = arb(10)**6
worstB = arb(0); worstcell=None; n=0
x = 400.0
while x < 100000.0:
    x2 = min(100000.0, x*1.02)
    v = LHS_cell(x, x2)
    assert v.upper() < C2, (x, x2, v.upper())
    if v.upper() > worstB: worstB = v.upper(); worstcell=(x,x2)
    n += 1
    x = x2
print(f"region B: {n} cells, worst LHS_cell upper = {worstB.str(10)} at {worstcell}, margin {float(1e6/worstB):.3f}x")
print(f"elapsed {time.monotonic()-t0:.1f}s")
open('progress.log','a').write(f"step: reran region B sweep, {n} cells, worst {worstB.str(6)}\n")
