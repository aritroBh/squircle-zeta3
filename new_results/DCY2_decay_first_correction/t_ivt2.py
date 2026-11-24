
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
N = arb(404, 4)
y = y_of(N); D = Delta_of(N)
eLam = (-kap_of(N)**2/4 + kap_of(N)*arb("1.117")*y.sqrt()).exp()
print("eLam:", eLam.str(8))
print("CH_R:", CH_R_of(N).str(8))
print("rDelta:", rDelta_of(N).str(8))
print("lamL:", lamL_of(N).str(8))
print("LogM(-D):", LogM_of(N,-D).str(8))
print("exp LogM(-D):", LogM_of(N,-D).exp().str(8))
print("CH_L/lamL:", ((PI/2).sqrt()*2**(arb(3)/4)/lamL_of(N)).str(8))
