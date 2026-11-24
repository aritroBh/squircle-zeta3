import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import flint
from flint import arb, ctx
from vlib import V_series
flint.ctx.cap = 64
ctx.dps = 110
t0 = time.time()
Vf, Vh = V_series(0, 12, h=arb(1)/32, T=4.6)
t1 = time.time()
print(f"V nlen=12 h=1/32: {t1-t0:.1f}s")
# extrapolate: nlen=41 at h=1/32 ~ (41/12)^1.5-ish x; h=1/64 ~ 4x nodes
