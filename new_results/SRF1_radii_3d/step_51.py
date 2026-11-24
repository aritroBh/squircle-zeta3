
import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
from flint import arb, arb_series, ctx
from vlib import E_integrand, ser, coeffs
ctx.dps = 60
X = ser(0, 12)
a = arb('0.3')
f = E_integrand(a, X, 1, 12)
print("prec:", f.prec)
cs = f.coeffs()
print("ncoeffs:", len(cs))
for i,c in enumerate(cs): print(i, c.str(10))
# check each factor's series length
la = a.log(); l1a = (1+a).log()
T1 = (X*(-1*l1a)).exp()/(1+a)
T2 = (1 + a*a*(X*(-2*la)).exp()).sqrt()
print("T1 prec", T1.prec, "T2 prec", T2.prec)
P = (X*(1*la)).exp()
print("P prec", P.prec)
print("T1 coeffs", len(T1.coeffs()), "T2", len(T2.coeffs()), "P", len(P.coeffs()))
