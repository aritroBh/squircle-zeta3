# Is sqrt2*L(chi_8,5) in the Q-span of the 65 sqrt2-free weight-5 atoms (products + Li5 at 8 points)?  (weight 3: 64 L(chi8,3) = sqrt2 (24 Li3(b) - 24 Li3(-b) - 4 L^3 + 3 pi^2 L))
import sys
sys.dont_write_bytecode = True
from id_lib import *
setprec(3100)
C = Consts()
D = 3000
tot = arb(0)
for r, ch in ((1, 1), (3, -1), (5, -1), (7, 1)):
    tot += ch * arb(5).zeta(arb(r) / 8)
L85 = tot / arb(8) ** 5
print("flint ctx.prec =", ctx.prec, "D =", D, "L(chi8,5) =", L85.str(25, radius=False), flush=True)
atoms = atoms5(C, with_sqrt2=False)
best, _ = lll_relation([L85 * C.sq2] + [v for _, v in atoms], D)
coeffs, last, mx = best
print("65 atoms, floor ~", f"{10**(D/66):.2e}", "maxcoef", f"{mx:.3e}", "last", last)
print({n: c for n, c in zip(["target"] + [n for n, _ in atoms], coeffs) if c != 0})
