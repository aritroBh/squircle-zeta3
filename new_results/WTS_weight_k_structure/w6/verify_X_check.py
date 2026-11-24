"""verify step: the 'total-derivative' integrals quoted in identify_w6.md,
   X_key = N[key|0] - N[key|a] - N[key|-b] = int_b^a (1/v - 1/(v-a) - 1/(v+b)) log^i t log^j(1-t) dv
         = int_0^1 log^i(t) log^j(1-t) dt/(1-t)   (because dlog(1-t) = (1/(v-a)+1/(v+b)-1/v) dv, dt = -(1-t) dlog(1-t)... sign checked numerically)
and the claimed closed forms X_50 = -8 pi^6/63, X_41 = 12 zeta3^2 - 2 pi^6/105, X_32 = 12 zeta3^2 - 23 pi^6/1260.
Evaluated BOTH as the v-integral (own Gauss-Legendre, exact distances) and as the t-integral (mp.quad, independent code path)."""
import sys, os, math
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp
from fractions import Fraction as Fr
DPS = 60
VL.set_dps(DPS + 10)
print("mp.mp.dps =", mp.mp.dps)
fin = {}
for (i, j) in ((5, 0), (4, 1), (3, 2)):
    fin[(i, j)] = {"0": (Fr(1), Fr(0)), "a": (Fr(-1), Fr(0)), "-b": (Fr(-1), Fr(0))}
vals = VL.quad_pieces(fin, list(fin), DPS + 10)
z3 = mp.zeta(3); pi6 = mp.pi ** 6
claimed = {(5, 0): -8 * pi6 / 63, (4, 1): 12 * z3 ** 2 - 2 * pi6 / 105, (3, 2): 12 * z3 ** 2 - 23 * pi6 / 1260}
for key in fin:
    i, j = key
    tq = mp.quad(lambda t: mp.log(t) ** i * mp.log1p(-t) ** j / (1 - t), [0, 0.5, 1])
    print(f"key {key}: v-integral (own GL) = {mp.nstr(vals[key], 30)}; claimed closed form diff = {mp.nstr(vals[key] - claimed[key], 4)}; "
          f"t-integral (mp.quad) - v-integral = {mp.nstr(tq - vals[key], 4)}   (mp.mp.dps={mp.mp.dps})")
