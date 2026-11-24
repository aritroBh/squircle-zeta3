
from flint import ctx, arb, acb
ctx.dps = 40
from step_07_partb import *
N = arb(404, 4)
print("y =", y_of(N).str(8))
print("S =", S_of(N).str(8))
print("exp S =", S_of(N).exp().str(8))
print("H =", H_of(N).str(8))
ms = momentsum_of(N)
print("momentsum =", ms.str(8))
print("ModelRatio =", ModelRatio_of(N).str(8))
print("T1ratio =", T1ratio_of(N).str(8))
# digamma/polygamma intervals
s0a = OMEGA/y_of(N)
print("digamma s0/2 =", (s0a/2).digamma().str(8))
print("polygamma s0/2 k=1 =", (s0a/2).polygamma(1).str(8))
print("polygamma s0/2 k=4 =", (s0a/2).polygamma(4).str(8))
