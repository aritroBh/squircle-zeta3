
import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import E_series, V_series, tompf

ctx.dps = 80; mp.mp.dps = 80
t0 = time.time()
Ef, Eh = E_series(2, 0, 6, h=arb(1)/32, T=4.6)
E = [tompf(x, 90) for x in Ef]
t1 = time.time()
Vf, Vh = V_series(0, 5, h=arb(1)/32, T=4.6)
V = [tompf(x, 90) for x in Vf]
t2 = time.time()
print(f"E: {t1-t0:.1f}s, V: {t2-t1:.1f}s")
s = [mp.mpf(24)] + [48*(E[k-1] + (V[k-2] if k>=2 else 0)) for k in range(1,6)]
claim = {1:'-26.558558386584754433798785874062052986842743535978204080483978',
         2:'-4.6892155184605022499600569810270175081392278040492906513172203885',
         3:'48.287432486014447759215859342074924467116625727590718495681774536',
         4:'-78.06628953961050950490562401928703057980'}
for k in range(1,6):
    line = f"s{k} = {mp.nstr(s[k], 50)}"
    if k in claim: line += f"  diff={mp.nstr(mp.mpf(claim[k])-s[k],3)}"
    print(line)
print("s5 =", mp.nstr(s[5], 30))
