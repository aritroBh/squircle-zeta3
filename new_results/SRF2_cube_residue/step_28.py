import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import W_int, J_int, tomp
from vlib import V_series, tompf
ctx.dps = 60; mp.mp.dps = 60

# sanity at eps=0.1 vs raw
Vf, _ = V_series(0.1, 2, h=arb(1)/32, T=4.6)
raw01 = tompf(Vf[0], 90)
eps = arb('0.1')
W = W_int(eps); J = J_int(eps); q = 1+2*eps
print("raw V(0.1)  =", mp.nstr(raw01, 20))
print("W+Pc (0.1)  =", mp.nstr(tomp(W - q*J/(eps+1)), 20), "  J(0.1) =", mp.nstr(tomp(J), 15))

# residue run
for e in [-0.9, -0.95, -0.99]:
    eps = arb(e)
    W = W_int(eps); J = J_int(eps); q = 1+2*eps
    Vv = W - q*J/(eps+1)
    print(f"eps={e}: W={mp.nstr(tomp(W),12)}  J={mp.nstr(tomp(J),12)}  V={mp.nstr(tomp(Vv),12)}  (1+eps)V={mp.nstr(tomp((1+eps)*Vv),15)}")

R_exact = -1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("R_exact =", mp.nstr(R_exact, 30))
