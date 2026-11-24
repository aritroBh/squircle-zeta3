import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import V_cont2, J2_int, W2_int, tomp
from vlib import V_series, tompf
ctx.dps = 60; mp.mp.dps = 60

# sanity at eps=0.1 vs raw
Vf, _ = V_series(0.1, 2, h=arb(1)/32, T=4.6)
raw01 = tompf(Vf[0], 90)
print("raw V(0.1)   =", mp.nstr(raw01, 20))
print("V_cont2(0.1) =", mp.nstr(tomp(V_cont2(arb('0.1'))), 20))

print("\nresidue run (improved):")
for e in ['-0.9', '-0.99', '-0.999']:
    eps = arb(e)
    Vv = V_cont2(eps)
    print(f"eps={e}: V={mp.nstr(tomp(Vv),15)}  (1+eps)V={mp.nstr(tomp((1+eps)*Vv),18)}")
print("\nR candidates: -1+I =", mp.nstr(-1+ (1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2),20),
      "  -2+I =", mp.nstr(-2+ (1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2),20))
