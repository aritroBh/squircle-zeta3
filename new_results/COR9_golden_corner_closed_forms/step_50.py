import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_lib import QuadField, QuadPoint, li4_coords, omega_class

QF5 = QuadField(5)
# Li4 args of the reduced basis as NF elements (a + b sqrt5)/e -> pair (Fr(a,e), Fr(b,e))
ARGS5 = {
 "Li4(-sqrt5 phi)": (Fr(-1,2), Fr(-1,2)),      # (-1-sqrt5)/2 = -phi... wait -sqrt5*phi = (-5-sqrt5)/2
 "Li4(1/(2phi^2))": (Fr(3,4), Fr(-1,4)),
 "Li4(1/phi^2)": (Fr(3,2), Fr(-1,2)),
 "Li4(1/2)": (Fr(1,2), Fr(0)),
 "Li4(2/phi^2)": (Fr(3), Fr(-1)),
 "Li4(phi/2)": (Fr(1,4), Fr(1,4)),
 "Li4(1/phi^6)": (Fr(9), Fr(-4)),
 "Li4(-1/phi^3)": (Fr(2), Fr(-1)),
}
ARGS5["Li4(-sqrt5 phi)"] = (Fr(-5,2), Fr(-1,2))  # (-5-sqrt5)/2
import mpmath as mp
mp.mp.dps = 50
r5 = mp.sqrt(5)
for nm, x in ARGS5.items():
    v = mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5
    print(nm, "=", mp.nstr(v, 12))

def lie_check(tag, c0, coefs, expect):
    P = QuadPoint(QF5, tag, c0, primes=[2])
    comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
    for nm, (p, q) in coefs.items():
        v = li4_coords(P, ARGS5[nm])
        for j in range(3):
            comp0[j] += p*v[j]; comp1[j] += q*v[j]
    cls = omega_class(P)
    ok = tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1])
    print(f"{tag}: Li4 part comps ({tuple(comp0)}, {tuple(comp1)})")
    print(f"   exact class(Omega)   ({tuple(cls[0])}, {tuple(cls[1])})   MATCH: {ok}; COR7 printed: {expect}")

lie_check("1/phi", (Fr(-1,2), Fr(1,2)), {
    "Li4(1/(2phi^2))": (Fr(-5,2), Fr(11,6)), "Li4(1/phi^2)": (Fr(113,16), Fr(-187,96)),
    "Li4(1/2)": (Fr(3,10), Fr(19,30)), "Li4(2/phi^2)": (Fr(-7,2), Fr(3,2)),
    "Li4(phi/2)": (Fr(-2), Fr(-2,3)), "Li4(1/phi^6)": (Fr(-29,72), Fr(1,12)),
    "Li4(-1/phi^3)": (Fr(35,18), Fr(7,6))}, "((111/2,0,-12), (-51/2,0,12))")
lie_check("1/phi^3", (Fr(-2), Fr(1)), {
    "Li4(1/(2phi^2))": (Fr(1,6), Fr(1,2)), "Li4(1/phi^2)": (Fr(-163,48), Fr(13,16)),
    "Li4(1/2)": (Fr(29,30), Fr(-11,10)), "Li4(2/phi^2)": (Fr(11,2), Fr(1,2)),
    "Li4(phi/2)": (Fr(11,3), Fr(1)), "Li4(1/phi^6)": (Fr(-1,8), Fr(1,72)),
    "Li4(-1/phi^3)": (Fr(1), Fr(7,9))}, "((84,-21,6), (-18,3,0))")
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_46: golden Lie-class checks.\n")
