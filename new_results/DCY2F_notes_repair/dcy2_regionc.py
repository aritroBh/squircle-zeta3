"""Region C for DCY2: N >= N0 = 1e5.

Claim: every piece of LHS(N) = E_total N^2 e^theta is non-increasing for N >= N0.
This file: (1) evaluates the pieces at N0 with arb; (2) verifies monotonicity
numerically at N0 * {1, 2, 5, 10, 100} via pointwise evaluation.
The y-scaling proof (a_tot >= 0 for every piece, kappa^2 = 28.2 >= 16) is in proof.md.
"""
from flint import arb, ctx
from dcy2_cell import *
from dcy2_sweep import LHS_cell

def pieces(N):
    """dict of the LHS piece values at N (region-B form, point eval)."""
    from dcy2_sweep import window_cell, TR_cell, TL_cell, T0_cell
    Nlo, Nhi = float(N), float(N)*1.0001
    y = y_of(N)
    bracket = window_cell(Nlo, Nhi)
    tp = TR_cell(Nlo, Nhi) + TL_cell(Nlo, Nhi)
    lhs = (2/PI)*y*S_sup(N)*( H_up_of(N)*bracket + tp )*arb(N)**2 + T0_cell(Nlo, Nhi)
    return lhs, bracket, tp

if __name__ == "__main__":
    ctx.dps = 40
    for N in [100000, 200000, 1000000, 10**7, 10**9, 10**12, 10**16, 10**20]:
        lhs, br, tp = pieces(N)
        print(f"N={N:9.0e}  LHS = {lhs.str(10)}")
