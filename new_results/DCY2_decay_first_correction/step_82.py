import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
cellcode = '''"""dcy2_sweep.py — the full cell sweep for DCY2.

Region A [400, 3000]: interval-N cells (ratio 1.02), window pieces by interval
arb, tails by endpoint monotonicity.
Region B [3000, 1e5]: endpoint monotonicity cells (ratio ~1.3).
Region C [1e5, inf): analytic single-point regime (monotone-decrease lemma).

Every value is a certified arb ball; cell accepted iff upper(LHS_cell) < C2.
"""
from flint import arb, acb, ctx
from dcy2_cell import *

C2 = arb(10)**6

def Gsplit_iv(jj, gam, c2m, D, ReL1, N):
    """Gsplit with interval quantities; take the max of both branches."""
    base_wide = None
    # branch 1: gam*D <= c2m/2 (check with interval compare)
    in_branch1 = gam*D <= c2m/2
    # branch 2 pieces
    T = min(D, c2m/(2*gam))
    # ReL1 here is ReL1max (function); Mtilde uses interval N via dcy2_cell
    if in_branch1:
        return Mtilde_of(N, jj, c2m/2)
    base = (gam*T**3).exp()*Mtilde_of(N, jj)
    excess = gam*D - c2m
    e2 = (ReL1*D).exp()
    if excess < 0:
        edge = EdgeInt_of(N, jj, -excess, T)*e2
    else:
        edge = e2*(excess*D**2).exp()*D**(jj+1)/(jj+1)
    return base + edge

def E_mom_iv(N):
    """E_mom with interval N (dcy2_cell functions accept arb intervals)."""
    return E_mom_of(N)

def E_tailmom_iv(N):
    return E_tailmom_of(N)

def tailpieces_cell_A(Nlo, Nhi):
    """region A: rates monotone increasing (min at Nlo), e^{-kap^2/4} max at Nlo."""
    eLam = (-kap_of(Nlo)**2/4 + kap_of(Nhi)*arb("1.117")*y_of(Nlo).sqrt()).exp()
    TR = eLam*CH_R_of(Nlo)/rDelta_of(Nlo)
    # left tail: Lambda(-Delta) <= |ReL1|Delta - c2 Delta^2 + Re(-l3 D^3 + l4 D^4 - l5 D^5 + l6 D^6) + B7tail
    y = y_of(Nlo); D = Delta_of(Nlo)
    lam_left = (ReL1max_of(Nlo)*D - 0*y, )  # placeholder
    return TR

def LHS_cell_A(Nlo, Nhi):
    """interval cell bound for region A."""
    N = arb((Nlo+Nhi)/2, (Nhi-Nlo)/2)
    bracket = E_mom_iv(N) + E_tailmom_iv(N)
    # tails: endpoint monotone
    eLam = (-kap_of(Nlo)**2/4 + kap_of(Nhi)*arb("1.117")*y_of(Nlo).sqrt()).exp()
    tailpieces = eLam*CH_R_of(Nlo)/rDelta_of(Nlo)
    # left tail pieces
    y = y_of(Nlo); D = Delta_of(Nlo)
    return (2/PI)*y*S_sup(Nlo)*( H_up_of(Nlo)*bracket + tailpieces )*arb(Nhi)**2

if __name__ == "__main__":
    ctx.dps = 40
    for lo, hi in [(400, 408), (500, 510), (1000, 1020), (2000, 2040), (2900, 2960)]:
        v = LHS_cell_A(lo, hi)
        print(f"[{lo},{hi}]  LHS_cell={v.str(10)}")
'''
open("dcy2_sweep.py","w").write(cellcode)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","dcy2_sweep.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
