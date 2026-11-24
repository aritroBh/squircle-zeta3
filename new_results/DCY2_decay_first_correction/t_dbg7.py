
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
for N in [10000, 10100, 10200, 100000, 101000]:
    N = arb(N)
    print(f"N={float(N):.0f}  LogM(-D) = {LogM_of(N,-Delta_of(N)).str(10)}")
Nlo, Nhi = 10000, 10200
A1, R_b = S.deriv_pieces(Nlo, Nhi)
print("A1 =", A1.str(8), " R_b =", R_b.str(8))
print("A0(mid) =", S.A0_of(arb((Nlo+Nhi)/2)).str(10))
