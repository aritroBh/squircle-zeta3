import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 1120
F4 = mp.mpf(open("hex_F4.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
def mpf(fr): return mp.mpf(fr.numerator)/fr.denominator
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2

# relation 1: basis Li4{1/3,2/3,3/4,1/2,1/9}
rel1 = {  # name -> Fr coeff
 "Li4(1/3)": Fr(103,8), "Li4(2/3)": Fr(-29,2), "Li4(3/4)": Fr(-21,4), "Li4(1/2)": Fr(-31,6), "Li4(1/9)": Fr(-45,64),
 "Li3(1/3)*log2": Fr(15,2), "Li3(-1/2)*log2": Fr(37,2), "Li3(-1/2)*log3": Fr(-5,4),
 "zeta3*log2": Fr(119,16), "zeta3*log3": Fr(-37,32),
 "Li2(1/3)^2": Fr(1), "Li2(1/3)*log2^2": Fr(-23,4), "Li2(1/3)*log2*log3": Fr(7,4), "Li2(1/3)*log3^2": Fr(1),
 "pi2*Li2(1/3)": Fr(-1,6), "pi2*log2^2": Fr(217,144), "pi2*log2*log3": Fr(-19,8), "pi2*log3^2": Fr(21,32),
 "pi4": Fr(83,540), "log2^4": Fr(-883,144), "log2^3*log3": Fr(269,24), "log2^2*log3^2": Fr(-57,8),
 "log2*log3^3": Fr(49,24), "log3^4": Fr(-55,192)}
# relation 2: basis Li4{-2,1/3,2/3,1/2,1/9}
rel2 = {
 "Li4(-2)": Fr(-21), "Li4(1/3)": Fr(-65,8), "Li4(2/3)": Fr(13,2), "Li4(1/2)": Fr(137,6), "Li4(1/9)": Fr(39,64),
 "Li3(1/3)*log2": Fr(15,2), "Li3(-1/2)*log2": Fr(37,2), "Li3(-1/2)*log3": Fr(-5,4),
 "zeta3*log2": Fr(119,16), "zeta3*log3": Fr(-37,32),
 "Li2(1/3)^2": Fr(1), "Li2(1/3)*log2^2": Fr(-23,4), "Li2(1/3)*log2*log3": Fr(7,4), "Li2(1/3)*log3^2": Fr(1),
 "pi2*Li2(1/3)": Fr(-1,6), "pi2*log2^2": Fr(-455,144), "pi2*log2*log3": Fr(9,8), "pi2*log3^2": Fr(-7,32),
 "pi4": Fr(-949,2160), "log2^4": Fr(-211,144), "log2^3*log3": Fr(101,24), "log2^2*log3^2": Fr(-15,8),
 "log2*log3^3": Fr(-35,24), "log3^4": Fr(113,192)}
vals = {
 "Li4(-2)": Li(4, mp.mpf(-2)), "Li4(1/3)": Li(4, mp.mpf(1)/3), "Li4(2/3)": Li(4, mp.mpf(2)/3),
 "Li4(3/4)": Li(4, mp.mpf(3)/4), "Li4(1/2)": Li(4, mp.mpf(1)/2), "Li4(1/9)": Li(4, mp.mpf(1)/9),
 "Li3(1/3)*log2": l313*l2, "Li3(-1/2)*log2": l31m*l2, "Li3(-1/2)*log3": l31m*l3,
 "zeta3*log2": z3*l2, "zeta3*log3": z3*l3,
 "Li2(1/3)^2": l213**2, "Li2(1/3)*log2^2": l213*l2**2, "Li2(1/3)*log2*log3": l213*l2*l3,
 "Li2(1/3)*log3^2": l213*l3**2, "pi2*Li2(1/3)": pi2*l213,
 "pi2*log2^2": pi2*l2**2, "pi2*log2*log3": pi2*l2*l3, "pi2*log3^2": pi2*l3**2, "pi4": mp.pi**4,
 "log2^4": l2**4, "log2^3*log3": l2**3*l3, "log2^2*log3^2": l2**2*l3**2, "log2*log3^3": l2*l3**3, "log3^4": l3**4}
for tag, rel in [("rel1", rel1), ("rel2", rel2)]:
    tot = mp.mpf(0)
    for nm, cf in rel.items():
        tot += mpf(cf)*vals[nm]
    print(tag, "|RHS - F4| =", mp.nstr(abs(tot - F4), 3))
