from flint import arb, ctx
ctx.dps = 45
PI = arb.pi()
Crest = arb(10_600_000)
delta = 2 - PI.sqrt()
def hN(N):
    N = arb(N)
    return N**2 * (-delta*N.sqrt()).exp()
# N3: Crest * h(N) <= 1
lo, hi = arb(100), arb(10**7)
while (hi-lo) > 1:
    mid = (lo+hi)/2
    if Crest*hN(mid) <= 1: hi = mid
    else: lo = mid
N3 = hi.ceil()
print("N3 =", N3.str(8))
# verify the inequality direction just below/above N3
print("Crest*h(N3-1) =", (Crest*hN(N3-1)).str(8), " (>1 expected)")
print("Crest*h(N3)   =", (Crest*hN(N3)).str(8), " (<=1 expected)")
# max of h on [100, N3]: h'(N)=0 at sqrt N = 4/delta -> N = 16/delta^2
Nc = 16/delta**2
print("16/delta^2 =", Nc.str(10))
# h decreasing for sqrt N >= 4/delta i.e. N >= 308.6..; check max on [100,309.5] at 309
import math
best = arb(0); bn=None
for k in range(100, 310):
    v = hN(k)
    if v > best: best = v; bn = k
print("max h on integers [100,309] =", best.str(10), "at N=", bn)
# and h decreasing beyond: h(309) vs h(310)
print("h(309)=", hN(309).str(10), " h(310)=", hN(310).str(10))
C3 = arb(10)**6 + Crest*hN(309)
print("C3 =", C3.str(12))
print("C3/1e10 =", (C3/arb(10)**10).str(8))
open('progress.log','a').write(f"step: rest combination N3={N3.str(6)}, C3={C3.str(8)}\n")
