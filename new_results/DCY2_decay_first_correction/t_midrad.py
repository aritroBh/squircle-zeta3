
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
for lo, hi in [(400, 408), (2000, 2040)]:
    N = arb((lo+hi)/2, (hi-lo)/2)
    em = S.E_mom_iv(N)
    et = S.E_tailmom_iv(N)
    print(f"[{lo},{hi}]")
    print("  E_mom  mid/rad:", em.mid(), em.rad())
    print("  E_tail mid/rad:", et.mid(), et.rad())
    v = S.LHS_cell_A(lo, hi)
    print("  LHS mid/rad:", v.mid(), v.rad())
