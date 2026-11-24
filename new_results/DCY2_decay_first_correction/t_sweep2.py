
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
