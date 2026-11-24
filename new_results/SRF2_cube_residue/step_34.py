import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import V_cont2, W2_int, J2_int, tomp
ctx.dps = 60; mp.mp.dps = 60

# V at safe distances from -1 (mass-hiding negligible for eps >= -0.9)
vals = {}
for e in ['-0.5', '-0.6', '-0.7', '-0.8', '-0.9']:
    eps = arb(e)
    Vv = V_cont2(eps)
    vals[float(e)] = tomp(Vv)
    print(f"V({e}) = {mp.nstr(vals[float(e)], 20)}   (1+eps)V = {mp.nstr((1+float(e))*vals[float(e)], 18)}")
