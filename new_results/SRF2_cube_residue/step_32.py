import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import W2_int, J2_int, tomp
ctx.dps = 60; mp.mp.dps = 60

for e in ['-0.9', '-0.99', '-0.999']:
    eps = arb(e)
    q = 1+2*eps
    W = W2_int(eps); J = J2_int(eps)
    Pc = -q/(eps+1)*J
    print(f"eps={e}: W2={mp.nstr(tomp(W),15)}  (1+eps)W2={mp.nstr(tomp((1+eps)*W),12)}  J2={mp.nstr(tomp(J),15)}")
    print(f"          -qJ={mp.nstr(tomp(-q*J),15)}  Pc={mp.nstr(tomp(Pc),15)}")
