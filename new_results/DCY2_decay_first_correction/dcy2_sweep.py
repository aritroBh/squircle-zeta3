"""dcy2_sweep.py — full cell sweep for DCY2 (region B: [400, 1e5]).

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
        # base: min( e^{gam T^3} Mtilde(c2m), widened Mtilde(c2m/2) ) on [0,T];
        # edge [T,Delta]: damped incomplete-Gamma integral (excess<0) or
        # polynomial moment with the excess factor (excess>0).
        gamD_lo = gam_lo*D_lo
        T_lo = min(D_lo, c_lo/(2*gam_lo))
        gamT3 = min(gam_lo*D_lo**3, c_lo**3/(8*gam_lo**2))   # R5: gamma T^3 <= max of both forms
        base = min(gamT3.exp()*Mtilde_of(Nhi, jj),
                   Mtilde_of(Nhi, jj, c_hi/2))
        if gamD_lo <= c_lo/2:
            return base
        excess_hi = gamD_lo - c_hi
        e2 = (ReL1max_of(Nlo)*D_lo).exp()
        if excess_hi < 0:
            rate_lo = c_hi - gamD_lo
            edge = EdgeInt_of(Nlo, jj, rate_lo, T_lo)*e2
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
        # R5: TP_p decreasing for p <= 2 (erfc/Q_1 form) -> Nlo; increasing for p >= 3 -> Nhi
        Et += cf*TailMom_of(Nlo if p <= 2 else Nhi, p)
    return E + 2*Et   # two-sided window tail (R4)

# ---------------- tail pieces ----------------
def TR_cell(Nlo, Nhi):
    """right tail T_R bound on the cell."""
    eLam = (-kap_of(Nlo)**2/4 + kap_of(Nhi)*arb("1.117")*y_of(Nlo).sqrt()).exp()
    rD = min(rDelta_of(Nlo), rDelta_of(Nhi))
    return eLam*CH_R_of(Nlo)/rD

def lamL_cell(Nlo, Nhi):
    return min(lamL_of(Nlo), lamL_of(Nhi))

def dN_of(N):
    N = arb(N); kap = kap_of(N)
    return 2/(N*kap**2) + 1/(4*N)

def LamPrime(N):
    """d Lambda(-Delta(N))/dN, algebraically cancelled form (tight for intervals).
       Lambda = N*A0 - ln u + bR*Delta, A0 = (1/2) ln(u^2 Q(0)/Q),
       Lambda' = A0 + N*A0' + (-ln u)' + bR Delta dN,
       A0' = u_N/u - (1/2)(Q_N/Q - Q0_N/Q0)
           = [sqrt2 y (u+sqrt2 y) u_N/u - (sqrt2 u+2y) y_N]/(2Q) + (sqrt2+2y) y_N/(2 Q0)."""
    N = arb(N)
    y = y_of(N); D = Delta_of(N)
    t = -D
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    dN = dN_of(N)
    yN = -y/(2*N)
    tN = t*dN
    uN = yN*t + y*tN
    A0 = (u*u*Q0/Q).log()/2
    A0p = (SQ2*y*(u+SQ2*y)*uN/u - (SQ2*u+2*y)*yN)/(2*Q) + (SQ2+2*y)*yN/(2*Q0)
    D1p = -uN/u + bR*D*dN
    return A0 + N*A0p + D1p

def TL_cell(Nlo, Nhi):
    """left tail T_L1 bound on the cell: LM <= LogM(mid) + V*hw, V >= |Lambda'| on cell."""
    Nmid = arb((Nlo+Nhi)/2)
    hw = arb(Nhi-Nlo)/2
    # V: evaluate LamPrime at endpoints and mid (exact points), take max |.|;
    # the interval evaluation of LamPrime on the whole cell certifies coverage.
    Vpts = max(abs(LamPrime(Nlo)), abs(LamPrime(Nmid)), abs(LamPrime(Nhi)))
    Vcell = abs(LamPrime(arb(Nmid, hw)))
    V = max(Vpts, Vcell)
    LM = LogM_of(Nmid, -Delta_of(Nmid)) + V*hw
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
