import sys, time
sys.path.insert(0, "<project>/research/08_theory/DCY2F_notes_repair")
from flint import arb, acb, ctx
ctx.dps = 35
from step_06_finite import main_arb, X_of, residual_check
for N in [100, 250, 400]:
    t0=time.monotonic()
    M = main_arb(N); dt=time.monotonic()-t0
    print(f"N={N}: main_arb {dt:.2f}s  M={M.str(10)}")
# sanity: residual_check at N=381 gives |M+X| ~ 1.62 scaled
M, res, bound = residual_check(381)
print("residual_check(381): res/bound =", (res/bound).str(8), " -> scaled const =", (res/bound*arb(10)**6).str(8))
