
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *

def LHS_iv(Nlo, Nhi):
    N = arb((Nlo+Nhi)/2, (Nhi-Nlo)/2)
    y = y_of(N)
    bracket = E_mom_of(N) + E_tailmom_of(N)
    D = Delta_of(N)
    eLam = (-kap_of(N)**2/4 + kap_of(N)*arb("1.117")*y.sqrt()).exp()
    tailpieces = eLam*CH_R_of(N)/rDelta_of(N) \
        + LogM_of(N,-D).exp()*(PI/2).sqrt()*2**(arb(3)/4)/lamL_of(N)
    return (2/PI)*y*S_sup(N)*(H_up_of(N)*bracket + tailpieces)*N**2 + T_0_of(N)

for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
    v = LHS_iv(lo, hi)
    print(f"[{lo},{hi}]  mid={v.mid()}  rad={v.rad()}")
