
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
N = arb(100000)
y = y_of(N); D = Delta_of(N)
print("y =", y.str(12), " Delta =", D.str(12), " kap =", kap_of(N).str(8))
print("LogM(-D) =", LogM_of(N,-D).str(12))
print("lamL_of =", lamL_of(N).str(12))
print("LogM(+D) =", LogM_of(N,D).str(12))
uD = 1-y*D
QD = uD**2 + SQ2*y*uD + y**2
P = -2*bR*uD**3 - (PI+2)*y*uD**2 + (SQ2*(N-2)*y**2 - 2*bR*y**3)*uD + 2*(N-1)*y**3
print("P(uD) =", P.str(12), "  pi y/2 =", (PI*y/2).str(12))
print("uD =", uD.str(10), " QD =", QD.str(10))
