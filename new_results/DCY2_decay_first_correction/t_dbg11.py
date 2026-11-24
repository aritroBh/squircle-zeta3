
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
lo, hi = 7799.841106705034, 7955.837928839134
br = S.window_cell(lo, hi)
tp = S.TR_cell(lo,hi) + S.TL_cell(lo,hi)
print("window bracket:", br.str(10))
print("pointwise window:", (E_mom_of(arb((lo+hi)/2)) + E_tailmom_of(arb((lo+hi)/2))).str(10))
print("tailpieces:", tp.str(10))
y = y_of(lo)
print("prefactor*(Hup*br+tp)*Nhi^2:", ((2/PI)*y*S_sup(lo)*(H_up_of(lo)*br+tp)*arb(hi)**2).str(10))
al = alpha_of(lo); B3 = B3_of(lo)
c_lo, c_hi = c2min_of(lo), c2min_of(hi)
D_lo, D_hi = Delta_of(lo), Delta_of(hi)
print("al*D_lo/(c_lo/2) =", (al*D_lo/(c_lo/2)).str(8))
print("B3*D_lo/(c_lo/2) =", (B3*D_lo/(c_lo/2)).str(8))
# alpha4 piece with the new base
print("Mtilde(Nhi,12,c_hi/2) =", Mtilde_of(hi,12,c_hi/2).str(8))
gamD_lo = al*D_lo
excess_hi = gamD_lo - c_hi
print("excess_hi =", excess_hi.str(8))
T_lo = min(D_lo, c_lo/(2*al))
e2 = (ReL1max_of(lo)*D_lo).exp()
edge = e2*(excess_hi*D_hi**2).exp()*D_hi**(12+1)/(12+1)
print("alpha4 edge =", edge.str(8), "  alpha4/24*edge =", (al**4/24*edge).str(8))
print("alpha4/24*Mtilde =", (al**4/24*Mtilde_of(hi,12,c_hi/2)).str(8))
# E_a piece
print("E_a =", (M4p_of(lo)/24*Mtilde_of(hi,4,c_hi/2)).str(8))
gamD_loB = B3*D_lo
excessB = gamD_loB - c_hi
T_loB = min(D_lo, c_lo/(2*B3))
edgeB = e2*(excessB*D_hi**2).exp()*D_hi**(4+1)/(4+1) if excessB>0 else D_hi**5*(-(c_hi-gamD_loB)*T_loB**2).exp()*e2/5
print("B3 edge =", edgeB.str(8), "  E_a edge =", (M4p_of(lo)/24*edgeB).str(8))
