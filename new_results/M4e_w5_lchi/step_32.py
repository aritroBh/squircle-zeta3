import time, math, json
import mpmath as mp
from flint import acb, arb, ctx

DIG = 1650
t0 = time.monotonic()
mp.mp.dps = DIG + 30
h = mp.mpf(8)**-5 * (mp.zeta(5, mp.mpf(1)/8) + mp.zeta(5, mp.mpf(7)/8)
                     - mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8))
print(f"mpmath hurwitz: {time.monotonic()-t0:.1f}s")
Lstr = mp.nstr(h, DIG + 10)
print(Lstr[:80], "...")

# flint cross-check at ~1670 digits
ctx.prec = math.ceil((DIG + 40) * math.log2(10))
t0 = time.monotonic()
z = (acb.zeta(acb(5), acb(1)/8) + acb.zeta(acb(5), acb(7)/8)
     - acb.zeta(acb(5), acb(3)/8) - acb.zeta(acb(5), acb(5)/8)) / acb(8)**5
print(f"flint hurwitz: {time.monotonic()-t0:.1f}s; radius={z.real.rad()}")
mid = z.real.mid()
import re
fstr = mid.str(DIG + 5, radius=False)
print(fstr[:80], "...")
d = abs(mp.mpf(fstr) - h)
print("mpmath-flint diff:", mp.nstr(d, 6))
json.dump(dict(digits=DIG, value=mp.nstr(h, DIG+5), flint_mid=fstr, flint_radius=str(z.real.rad())),
          open('<project>/research/08_theory/M4e_w5_lchi/Lchi8_5.json','w'), indent=1)
print("saved Lchi8_5.json")
