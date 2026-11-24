import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import V_series, tompf
flint.ctx.cap = 64
ctx.dps = 60; mp.mp.dps = 60
# Confirm JSON V coefficients independently (order 8, h=1/32)
Vf, Vh = V_series(0, 8, h=arb(1)/32, T=4.6)
V8 = [tompf(x, 70) for x in Vf]
V8h = [tompf(x, 70) for x in Vh]
for k in range(8):
    print(k, mp.nstr(V8[k], 25), " hdiff:", mp.nstr(abs(V8[k]-V8h[k]), 3))
