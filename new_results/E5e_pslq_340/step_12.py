import json
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 340
C = json.load(open("consts360.json"))
R = json.load(open("relation.json"))
rel = R["rel"]; q = [Fraction(s) for s in R["q"]]
labels = C["labels"]
strs = [C["vals"][l] for l in labels]
v = [mp.mpf(s) for s in strs]
X = mp.mpf(C["vals"]["X"])

closed = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*vv for x, vv in zip(q, v))
print("X (data)       =", mp.nstr(X, 60))
print("closed form    =", mp.nstr(closed, 60))
print("difference     =", mp.nstr(X-closed, 6))

d = json.load(open('../E5d_nonlocal_pslq/b0_p5_convB_350.json'))
P5 = mp.mpf(d['e1']['P5'])
print("\nP5 (data)      =", mp.nstr(P5, 50))
print("sqrt5*closed   =", mp.nstr(mp.sqrt(5)*closed, 50))
print("difference     =", mp.nstr(P5 - mp.sqrt(5)*closed, 6))

# --- local law comparison ---
sq5 = mp.sqrt(5)
z5c  = q[labels.index("z5")]      # coord of zeta(5) in X = P/√5
z23c = q[labels.index("z2z3")]    # coord of zeta(2)zeta(3)
print("\ncoord[z5]  (X-space) =", z5c, "=", float(z5c))
print("coord[z2z3](X-space) =", z23c, "=", float(z23c))
loc5  = Fraction(-2079,250); loc23 = Fraction(-477,250)   # X-space local coords
N5x  = z5c - loc5;  N23x = z23c - loc23
print("\nX-space: N5  = coord - (-2079/250) =", N5x, "=", float(N5x))
print("X-space: N23 = coord - (-477/250)  =", N23x, "=", float(N23x))
# P-space (coords carry √5):  N^(5) = coord_P - (-2079√5/250), etc.
N5P  = N5x  # rational multiplier of sqrt5
print("\nP-space: N^(5)   = (%s)·√5 = %s√5  ≈ %s" % (N5P, N5P, mp.nstr(mp.mpf(N5P.numerator)/N5P.denominator*sq5, 20)))
print("P-space: N^(2,3) = (%s)·√5 = %s√5  ≈ %s" % (N23x, N23x, mp.nstr(mp.mpf(N23x.numerator)/N23x.denominator*sq5, 20)))
print("\nsimplified: N^(5) = -711√5/400 ?", Fraction(-711,400) == N5x)
print("simplified: N^(2,3) = -1239√5/250 ?", Fraction(-1239,250) == N23x)
