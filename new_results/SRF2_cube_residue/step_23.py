import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import V_series, tompf
ctx.dps = 80; mp.mp.dps = 80
import json
d = json.load(open("<project>/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json"))
V = [mp.mpf(x) for x in d["V"]]

# 1) raw V(0.1) from vlib series
Vf, _ = V_series(0.1, 2, h=arb(1)/32, T=4.6)
print("raw V(0.1)   =", mp.nstr(tompf(Vf[0],90), 20))

# 2) W + Pc at 0.1
eps = arb('0.1')
W = W_int(eps); J = J_int(eps); q = 1+2*eps
print("W+Pc at 0.1  =", mp.nstr(tomp(W - q/(eps+1)*J), 20))

# 3) V(-0.5): Taylor series sum (radius 1, converges at -0.5) vs W+Pc
tay = sum(V[k]*mp.mpf('-0.5')**k for k in range(31))
print("Taylor V(-0.5)=", mp.nstr(tay, 20))
eps = arb('-0.5')
W = W_int(eps); J = J_int(eps); q = 1+2*eps
print("W+Pc at -0.5 =", mp.nstr(tomp(W - q/(eps+1)*J), 20))
