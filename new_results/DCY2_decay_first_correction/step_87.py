import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
sweep = '''"""dcy2_sweep.py — full cell sweep for DCY2 (region B: [400, 1e5]).

Design (D13-style):
- narrow cells (ratio 1.02); every factor assigned to a cell endpoint by its
  monotonicity class (table in proof.md);
- tails by endpoint monotonicity; the single problematic factor LogM(-Delta)
  is enclosed by the mean-value form  LogM(mid) + V * halfwidth  with the
  signed cancellation  ln u - (1/2) ln(Q/Q0)  kept exact in the derivative;
- region C [1e5, inf): analytic regime (y-scaling, kappa^2 = 28.2 >= 16).
"""
from flint import arb, acb, ctx
from dcy2_cell import *

C2 = arb(10)**6

# ---------------- window pieces with endpoint assignment ----------------
def window_cell(Nlo, Nhi):
    """E_mom + E_tailmom for cell [Nlo, Nhi]. Every factor at its monotone endpoint:
       decreasing -> Nlo, increasing -> Nhi."""
    y_lo, y_hi = y_of(Nlo), y_of(Nhi)
    D_lo, D_hi = Delta_of(Nlo), Delta_of(Nhi)
    c_lo, c_hi = c2min_of(Nlo), c2min_of(Nhi)
    # alpha, beta2, beta3, L_m, H_j: decreasing -> Nlo
    al = alpha_of(Nlo); be2 = beta2_of(Nlo); be3 = beta3_of(Nlo)
    L3,L4,L5 = L_of(Nlo,3),L_of(Nlo,4),L_of(Nlo,5)
    B3,B6 = B3_of(Nlo), B6_of(Nlo)
    H1 = arb("1.5")*y_lo*(1+2*y_lo); H2 = arb("1.875")*y_lo**2*(1+2*y_lo); H3 = arb("2.1875")*y_lo**3*(1+2*y_lo)
    M4 = M4p_of(Nlo)
    Hsup = 1 + H1*D_hi + H2*D_hi**2 + H3*D_hi**3
    # Gsplit pieces: branch check with the (decreasing) ratio gam*Delta/(c2min/2):
    # ratio decreasing in N (kappa*y^{3/2}-family) -> if holds at Nlo it holds on the cell
    def GsplitB(jj, gam_lo):
        # gam = B3 or alpha evaluated at Nlo (both decreasing)
        gamD_lo = gam_lo*D_lo
        if gamD_lo <= c_lo/2:
            # branch 1 on the whole cell: widened moment, c2min(Nhi)/2
            return Mtilde_of(Nhi, jj, c_hi/2)
        # branch 2: base + edge, endpoint rules
        T_lo = min(D_lo, c_lo/(2*gam_lo))
        gamT3 = gam_lo*T_lo**3
        base = gamT3.exp()*Mtilde_of(Nhi, jj)
        excess_hi = gamD_lo - c_hi
        e2 = (ReL1max_of(Nlo)*D_lo).exp()
        if excess_hi < 0:
            rate_lo = c_hi - gamD_lo
            edge = D_hi**(jj+1)*(-rate_lo*T_lo**2).exp()*e2/(jj+1)
        else:
            edge = e2*(excess_hi*D_hi**2).exp()*D_hi**(jj+1)/(jj+1)
        return base + edge
    E = (M4/24)*GsplitB(4, B3)
    E += Hsup*( be2*Mtilde_of(Nhi,8) + be3*Mtilde_of(Nhi,10) + al**4/24*GsplitB(12,al)
                + B6*(B6*D_lo**6).exp()*GsplitB(6,al) )
    E_drop = {6: H1*L5+H2*L4+H3*L3, 7: H2*L5+H3*L4,
              8: H1*L3*L4+H2*L3**2/2+H3*L5, 9: H2*L3*L4+H3*L3**2/2,
              10: H1*L3**3/6, 11: H2*L3**3/6, 12: H3*L3**3/6}
    E += sum(cf*Mtilde_of(Nhi,p) for p,cf in E_drop.items())
    # E_tailmom
    pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
          6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
    Et = arb(0)
    for p, cf in pp.items():
        Et += cf*TailMom_of(Nhi, p)
    return E + Et

# ---------------- tail pieces ----------------
def TR_cell(Nlo, Nhi):
    """right tail T_R bound on the cell."""
    eLam = (-kap_of(Nlo)**2/4 + kap_of(Nhi)*arb("1.117")*y_of(Nlo).sqrt()).exp()
    rD = min(rDelta_of(Nlo), rDelta_of(Nhi))
    return eLam*CH_R_of(Nlo)/rD

def lamL_cell(Nlo, Nhi):
    return min(lamL_of(Nlo), lamL_of(Nhi))

def LogM_mD_mv(Nmid, V, hw):
    """LogM(-Delta) at Nmid plus mean-value allowance V*hw."""
    return LogM_of(Nmid, -Delta_of(Nmid)) + V*hw

def dLam_dN(N):
    """V >= |d/dN LogM(-Delta(N))| on cell containing N (interval eval OK).
       d/dN LogM(t;N) = ln u - (1/2) ln(Q/Q0) + (N-1) u_N/u - (N/2)(Q_N/Q - Q0_N/Q0) - bR t dN,
       t = -Delta, dN = d ln Delta/dN = 2/(N kap^2) + 1/(4N)."""
    N = arb(N)
    y = y_of(N); kap = kap_of(N); D = Delta_of(N)
    t = -D
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    dN = 2/(N*kap**2) + 1/(4*N)
    yN = -y/(2*N)
    tN = t*dN
    uN = yN*t + y*tN
    QN = (2*u+SQ2*y)*uN + (SQ2*u+2*y)*yN
    Q0N = (SQ2+2*y)*yN
    D0 = (u).log() - (Q/Q0).log()/2
    D1 = (N-1)*uN/u - (N/2)*(QN/Q - Q0N/Q0) - bR*t*dN
    return abs(D0+D1)

def TL_cell(Nlo, Nhi):
    """left tail T_L1 bound on the cell (mean-value LogM)."""
    Nmid = (Nlo+Nhi)/2
    hw = (Nhi-Nlo)/2
    V = dLam_dN(arb(Nmid, hw))
    LM = LogM_mD_mv(Nmid, V, hw)
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    return LM.exp()*CH_L/lamL_cell(Nlo, Nhi)

def T0_cell(Nlo, Nhi):
    """[0,1] piece: T_0(N) decreasing for N >= 400 -> at Nlo, with e^{theta} at Nhi."""
    N = arb(Nlo)
    th_hi = (PI*arb(Nhi)).sqrt()
    base = (PI/2).sqrt()*2**(arb(3)/4)*(2**(arb(3)/2)*(-arb("0.5348")*(N-arb(5)/2)).exp()
        + 2**(-N+arb(5)/2)*2**(arb(3)/2))
    return base*arb(Nhi)**2*th_hi.exp()

# ---------------- full LHS for a cell ----------------
def LHS_cell(Nlo, Nhi):
    bracket = window_cell(Nlo, Nhi)
    tailpieces = TR_cell(Nlo, Nhi) + TL_cell(Nlo, Nhi)
    y = y_of(Nlo)
    return (2/PI)*y*S_sup(Nlo)*( H_up_of(Nlo)*bracket + tailpieces )*arb(Nhi)**2 + T0_cell(Nlo, Nhi)

def sweep(lo, hi, ratio=1.02, verbose=True):
    """sweep cells [lo, hi] geometric in ratio; return max LHS_cell upper."""
    worst = arb(0)
    x = lo
    while x < hi:
        x2 = min(hi, x*ratio)
        v = LHS_cell(x, x2)
        worst = max(worst, v.upper())
        if verbose:
            print(f"  cell [{x:.1f}, {x2:.1f}]  LHS upper = {v.upper().str(10)}")
        x = x2
    return worst

if __name__ == "__main__":
    ctx.dps = 40
    for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
        v = LHS_cell(lo, hi)
        print(f"[{lo},{hi}]  LHS upper = {v.upper().str(10)}  (mid {v.mid()})")
'''
open("dcy2_sweep.py","w").write(sweep)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","dcy2_sweep.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
