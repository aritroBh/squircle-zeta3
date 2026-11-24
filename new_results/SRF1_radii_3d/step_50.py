
import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
from flint import arb, arb_series, ctx
from vlib import E_integrand, ser, coeffs
ctx.dps = 60
X = ser(0, 12)
a = arb('0.3')
f = E_integrand(a, X, 1, 12)
print("type:", type(f), "prec:", f.prec() if hasattr(f,'prec') else '?')
cs = f.coeffs()
print("ncoeffs:", len(cs))
for i,c in enumerate(cs): print(i, c.str(12))
