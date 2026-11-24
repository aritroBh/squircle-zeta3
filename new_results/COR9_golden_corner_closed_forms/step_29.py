import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from fractions import Fraction as Fr
from cor9_lib import QuadField, QuadPoint, li4_coords, omega_class, L7

# ---- exact Lie class of the Li4 part of the W4 closed form, corner (sqrt2-1)^2 ----
QF = QuadField(2)
P = QuadPoint(QF, "(sqrt2-1)^2", (Fr(3), Fr(-2)), primes=[2])
args = {"b": (Fr(-1), Fr(1)), "-b": (Fr(1), Fr(-1)), "-b^2": (Fr(-3), Fr(2)),
        "1/sqrt2": (Fr(0), Fr(1, 2)), "1/2": (Fr(1, 2), Fr(0)),
        "1-1/sqrt2": (Fr(1), Fr(-1, 2)), "2-sqrt2": (Fr(2), Fr(-1))}
coef = {"b": (Fr(-478), Fr(-153)), "-b": (Fr(-470), Fr(-153)), "-b^2": (Fr(83,2), Fr(59,4)),
        "1/sqrt2": (Fr(4), Fr(-4)), "1/2": (Fr(103,2), Fr(41,4)),
        "1-1/sqrt2": (Fr(16), Fr(14)), "2-sqrt2": (Fr(116), Fr(30))}
comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
for nm, x in args.items():
    v = li4_coords(P, x)
    p, q = coef[nm]
    for j in range(3):
        comp0[j] += p*v[j]; comp1[j] += q*v[j]
cls = omega_class(P)
print("Li4 part of numerical W4:  comp0:", tuple(comp0), " comp1:", tuple(comp1))
print("class(Omega((sqrt2-1)^2)): comp0:", tuple(cls[0]), " comp1:", tuple(cls[1]))
print("MATCH:", tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1]))
print("COR7 printed: (84,-100,128), (-4,-16,48)")

# ---- W2 closed form ----
mp.mp.dps = 1120
W2v = mp.mpf(open(WD+"b2sq_W2_hpl_1150.txt").read())
b = mp.sqrt(2)-1; c = b*b
Li = lambda s, x: mp.polylog(s, x)
W2cf = -Li(2, c) - Li(2, -c) + mp.pi**2/12
W2cf2 = -Li(2, c**2)/2 + mp.pi**2/12   # by duplication Li2(c)+Li2(-c) = Li2(c^2)/2
print("|W2 - (-Li2(c)-Li2(-c)+pi^2/12)| =", mp.nstr(abs(W2v-W2cf), 3))
print("|W2 - (-Li2(c^2)/2+pi^2/12)|    =", mp.nstr(abs(W2v-W2cf2), 3))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_29: b2sq Lie-class check + W2 closed form.\n")
