"""RST1 step_10: validate the assembled analytic bound against the exact rest."""
import mpmath as mp
from rst1_common import *
rows = load_mk()
print(" k      |rest|            bound(k)         C_eff=bound/target   ratio |rest|/bound")
for k in [99,120,150,200,250,300,350,400]:
    mm = Mmain(k); rest = abs(rows[k]-mm)
    b = rest_bound(k); t = target(k)
    print(f"{k:4d}  {mp.nstr(rest,6):>12s}   {mp.nstr(b,6):>12s}   {mp.nstr(b/t,8):>10s}   {mp.nstr(rest/b,8)}")
# constant C as max of bound/target over k>=99 (sample)
mx=0; mk=None
for k in range(99,401,2):
    r=rest_bound(k)/target(k)
    if r>mx: mx=r; mk=k
print("max bound/target over k in [99,400] step2:", mp.nstr(mx,8), "at k=",mk)
