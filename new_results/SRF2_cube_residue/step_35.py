import sys, json
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import V_cont2, tomp
ctx.dps = 60; mp.mp.dps = 60
d = json.load(open("<project>/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json"))
V = [mp.mpf(x) for x in d["V"]]

for e in ['-0.3', '-0.6', '-0.7', '-0.8', '-0.9']:
    x = mp.mpf(e)
    tay = sum(V[k]*x**k for k in range(31))
    print(f"Taylor V({e}) = {mp.nstr(tay, 18)}")
print()
for e in ['-0.3', '-0.6', '-0.7', '-0.8', '-0.9']:
    v = tomp(V_cont2(arb(e)))
    print(f"V_cont2({e}) = {mp.nstr(v, 18)}")
