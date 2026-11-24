
from flint import ctx
ctx.dps = 40
from dcy2_sweep import LHS_cell, TL_cell, TR_cell
for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
    v = LHS_cell(lo, hi)
    print(f"[{lo},{hi}]  LHS upper = {v.upper().str(10)}   TL = {TL_cell(lo,hi).str(8)}")
