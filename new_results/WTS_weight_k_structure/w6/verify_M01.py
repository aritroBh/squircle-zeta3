"""verify step: own definition-based M_0, M_1 (verify_lib.direct_M) vs the PROVED closed forms of the paper
(M_0 = log2 - sqrt2 L, M_1 of thm:M1), L = log(1+sqrt2)."""
import sys, os, math
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp
TGT = 60
Y = math.log(10) * (TGT + 12) + 25
wd = TGT + int(Y / 2.302585) + 30
VL.set_dps(wd)
m = [VL.direct_M(k, TGT, Y=Y) for k in (0, 1)]
VL.set_dps(TGT + 10)
print("mp.mp.dps =", mp.mp.dps, "(working precision during direct_M:", wd, ")")
lam = mp.log(2); L = mp.log(1 + mp.sqrt(2)); b = mp.sqrt(2) - 1
M0 = lam - mp.sqrt(2) * L
M1 = (19 * mp.sqrt(2) / 48 - mp.mpf(1) / 24) * mp.pi ** 2 + lam ** 2 / 2 + 2 * mp.sqrt(2) * lam * L - (1 + 5 * mp.sqrt(2) / 2) * L ** 2 - 5 * mp.sqrt(2) * mp.polylog(2, b)
print("M_0 own - closed form =", mp.nstr(m[0] - M0, 5))
print("M_1 own - closed form =", mp.nstr(m[1] - M1, 5))
