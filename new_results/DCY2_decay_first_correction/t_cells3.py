
from flint import ctx, arb
ctx.dps = 40
from dcy2_sweep import LHS_cell, TL_cell, LamPrime
from dcy2_cell import LogM_of, Delta_of
for N0 in [404, 635, 2020, 10100, 101000]:
    h = max(N0//1000, 1)
    fd = (LogM_of(N0+h,-Delta_of(N0+h)) - LogM_of(N0-h,-Delta_of(N0-h)))/(2*h)
    print(f"N={N0}: FD={fd.str(10)}  LamPrime={LamPrime(N0).str(10)}")
print()
for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
    v = LHS_cell(lo, hi)
    print(f"[{lo},{hi}]  LHS upper = {v.upper().str(10)}   TL = {TL_cell(lo,hi).str(8)}")
