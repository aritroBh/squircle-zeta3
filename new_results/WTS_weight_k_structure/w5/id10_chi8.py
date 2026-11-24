# Is L(chi_8,5) (even real character mod 8, weight 5) and L(chi_8,3) in the span of the level-8 real polylog atoms?
import sys
sys.dont_write_bytecode = True
from id_lib import *
setprec(1500)
C = Consts()
def Lchi8(s):
    tot = arb(0)
    for r, ch in ((1, 1), (3, -1), (5, -1), (7, 1)):
        tot += ch * arb(s).zeta(arb(r) / 8)
    return tot / arb(8) ** s
D = 1400
print("flint ctx.prec =", ctx.prec, "D =", D)
L83 = Lchi8(3); L85 = Lchi8(5)
print("L(chi8,3) =", L83.str(30, radius=False), " L(chi8,5) =", L85.str(30, radius=False))
b = C.b; s2 = C.sq2; lam = C.lam; L = C.L; pi2 = C.pi**2
# weight 3 test: L(chi8,3) in span_Q(sqrt2){zeta3, Li3(b), Li3(-b), Li3(1/s2), lam^3,.., pi2 lam, pi2 L}
w3 = [("z3", C.z3), ("Li3b", C.Li(3, b)), ("Li3-b", C.Li(3, -b)), ("Li3(1/s2)", C.Li(3, 1/s2))] + [(f"lam^{i} L^{3-i}", lam**i * L**(3-i)) for i in range(4)] + [("pi2 lam", pi2*lam), ("pi2 L", pi2*L)]
w3 = w3 + [("sq2*" + n, v*s2) for n, v in w3]
best, _ = lll_relation([L83] + [v for _, v in w3], D)
print("L(chi8,3) vs weight-3 atoms (x{1,sqrt2}):", "maxcoef", best[2], {n: c for n, c in zip(["target"] + [n for n, _ in w3], best[0]) if c != 0})
# weight 5: Li5 at 8 points, z5, sqrt2 multiples ; also with products? only pure first
pts = {"1/2": arb(1)/2, "b": b, "-b": -b, "1/s2": 1/s2, "-b2": -b*b, "2-s2": 2-s2, "1-1/s2": 1-1/s2, "1/s2-1": 1/s2-1}
w5 = [("z5", C.z5)] + [(f"Li5({n})", z.polylog(5)) for n, z in pts.items()]
w5 = w5 + [("sq2*" + n, v*s2) for n, v in w5]
best, _ = lll_relation([L85] + [v for _, v in w5], D)
print("L(chi8,5) vs Li5 atoms (x{1,sqrt2}) [pure weight-5 gens only]:", "maxcoef", f"{best[2]:.3e}", {n: c for n, c in zip(["target"] + [n for n, _ in w5], best[0]) if c != 0})
