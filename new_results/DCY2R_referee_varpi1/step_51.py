import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
from dcy2_cell import *
import dcy2_sweep as S

# (1) TailMom monotonicity: p=0,1,2 decreasing? evaluate across the worst cell and one more
print("TailMom_of(N,p) across worst cell [3075.23, 3136.74]:")
for p in [0,1,2,3,9]:
    vlo = TailMom_of(3075.23, p); vhi = TailMom_of(3136.74, p)
    print(f"  p={p}: Nlo {vlo.str(8)}  Nhi {vhi.str(8)}  ratio hi/lo {float(vhi/vlo):.4f}")

# (2) quantify the E_tailmom(p<=2) share of LHS on the worst cell
def Et_share(N):
    y = y_of(N)
    H1 = arb("1.5")*y*(1+2*y); H2 = arb("1.875")*y**2*(1+2*y)
    return TailMom_of(N,0) + H1*TailMom_of(N,1) + H2*TailMom_of(N,2)
for N in [400.0, 3075.23, 3136.74, 100000.0]:
    e = Et_share(N)
    lhs = LHS_of(N)
    print(f"N={N:9.1f}: E_tailmom(p<=2) = {e.str(8)},  ~share of LHS = {float(e*1.0/lhs*1e6):.4f}e-6 * LHS (times (2/pi)y S H N^2 = {float((2/PI)*y_of(N)*S_sup(N)*H_up_of(N)*arb(N)**2):.3e})")
