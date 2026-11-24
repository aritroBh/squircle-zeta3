import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
import dcy2_cell as C
from dcy2_cell import *
import dcy2_sweep as S

# --- patch: Mtilde with provably-valid correction  2 efac tR^{j+1}/(j+1) ---
def Mtilde_valid(N, jj, c2use=None):
    c2m = C.c2min_of(N) if c2use is None else c2use
    ReL1 = C.ReL1max_of(N); tR = ReL1/(2*c2m)
    tot = arb(0)
    for i in range(jj+1):
        tot += C.binom(jj,i)*tR**(jj-i)*c2m**(-(arb(i)+1)/2)*((arb(i)+1)/2).gamma()
    efac = (ReL1**2/(4*c2m)).exp()
    return efac*tot + 2*efac*tR**(jj+1)/(jj+1)

# --- patch: TailMom x2 (provably valid two-sided bound) ---
def TailMom_valid(N, p):
    return 2*C.TailMom_of(N, p)

def window_cell_repaired(Nlo, Nhi):
    """mirror of dcy2_sweep.window_cell with repaired moments (pointwise at Nlo=Nhi=N)."""
    N = Nlo
    y = y_of(N); D = Delta_of(N)
    c2m = c2min_of(N)
    al = alpha_of(N); be2 = beta2_of(N); be3 = beta3_of(N)
    L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    B3,B6 = B3_of(N), B6_of(N)
    H1 = arb("1.5")*y*(1+2*y); H2 = arb("1.875")*y**2*(1+2*y); H3 = arb("2.1875")*y**3*(1+2*y)
    M4 = M4p_of(N)
    Hsup = 1 + H1*D + H2*D**2 + H3*D**3
    def GsplitR(jj, gam):
        if gam*D <= c2m/2:
            return Mtilde_valid(N, jj, c2m/2)
        T = min(D, c2m/(2*gam))
        base = (gam*T**3).exp()*Mtilde_valid(N, jj)
        excess = gam*D - c2m
        e2 = (ReL1max_of(N)*D).exp()
        if excess < 0:
            edge = EdgeInt_of(N, jj, -excess, T)*e2
        else:
            edge = e2*(excess*D**2).exp()*D**(jj+1)/(jj+1)
        return base + edge
    E = (M4/24)*GsplitR(4, B3)
    E += Hsup*( be2*Mtilde_valid(N,8) + be3*Mtilde_valid(N,10) + al**4/24*GsplitR(12,al)
                + B6*(B6*D**6).exp()*GsplitR(6,al) )
    E_drop = {6: H1*L5+H2*L4+H3*L3, 7: H2*L5+H3*L4,
              8: H1*L3*L4+H2*L3**2/2+H3*L5, 9: H2*L3*L4+H3*L3**2/2,
              10: H1*L3**3/6, 11: H2*L3**3/6, 12: H3*L3**3/6}
    E += sum(cf*Mtilde_valid(N,p) for p,cf in E_drop.items())
    pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
          6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
    Et = sum(cf*TailMom_valid(N,p) for p,cf in pp.items())
    return E + Et

def LHS_repaired(N):
    Nlo, Nhi = N, N*1.0001
    y = y_of(N)
    bracket = window_cell_repaired(N, N)
    tp = S.TR_cell(Nlo,Nhi) + S.TL_cell(Nlo,Nhi)
    return (2/PI)*y*S_sup(N)*( H_up_of(N)*bracket + tp )*arb(N)**2 + S.T0_cell(Nlo,Nhi)

for N in [400, 3075, 10000, 100000]:
    lr = LHS_repaired(float(N))
    lo = C.LHS_of(N)
    print(f"N={N:7d}: LHS_orig(point) = {lo.str(10)}   LHS_repaired = {lr.str(10)}   ratio {float(lr/lo):.4f}")
