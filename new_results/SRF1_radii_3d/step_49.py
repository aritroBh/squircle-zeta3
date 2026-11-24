
import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import E_series, tompf
ctx.dps = 110; mp.mp.dps = 110
E2f, _ = E_series(1, 0, 12, h=arb(1)/32, T=4.6)
for k, x in enumerate(E2f):
    print(k, x.mid().str(30, radius=False), "+-", x.rad().str(5, radius=False))
