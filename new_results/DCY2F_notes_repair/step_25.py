import sys, time
sys.path.insert(0, "<project>/research/08_theory/DCY2_decay_first_correction")
from flint import arb, acb, ctx
ctx.dps = 40
t0=time.monotonic()
from dcy2_sweep import LHS_cell
v = LHS_cell(3075.23, 3136.74)
print("current code  LHS_cell(3075.23,3136.74) =", v.upper().str(12))
# Also a couple of other cells near the top for context
for lo,hi in [(400,408),(1000,1020),(5000,5100),(10000,10200),(50000,51000),(99990,100000)]:
    vv = LHS_cell(lo,hi)
    print(f"  cell [{lo},{hi}] -> {vv.upper().str(10)}")
print("elapsed", time.monotonic()-t0)
