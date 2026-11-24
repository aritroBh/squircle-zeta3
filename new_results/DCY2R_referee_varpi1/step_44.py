import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
from dcy2_cell import LHS_of
t0 = time.monotonic()
prev = None; bad = 0; firstN=None
N = 100000.0
vals = []
while N < 1e7:
    v = LHS_of(N)
    if prev is not None and v.upper() > prev.upper():
        bad += 1
        if firstN is None: firstN = N
    prev = v
    N *= 1.005
print(f"dense check [1e5, 1e7] ratio 1.005: {int((7-5)/__import__('math').log10(1.005))} points, violations of non-increasing: {bad}, first at N={firstN}")
for N in [1e7, 1e9, 1e12, 1e16, 1e20]:
    print(f"  LHS_of(1e{int(__import__('math').log10(N))}) = {LHS_of(N).str(8)}")
print("elapsed", time.monotonic()-t0)
open('progress.log','a').write(f"step: region C dense monotonicity check, violations={bad}\n")
