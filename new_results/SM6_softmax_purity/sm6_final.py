# SM6: final closed forms, cross-checks and negative PSLQs for c6, c7.
import mpmath as mp, pickle
from fractions import Fraction as Fr
mp.mp.dps = 420
z = mp.zeta
CAND, cv, gams = pickle.load(open("sm6_w7_raw.pkl", "rb"))
g33, g15 = gams
X1, X2 = cv[1], cv[2]           # int (-log x)^3 log^2(1-x) log(1+x)/x ,  int (-log x)^3 log(1-x) log^2(1+x)/x
Xs = X1 + X2                    # = int (-log x)^3 log(1-x) log(1+x) log(1-x^2)/x dx
# weight 7 relation found by PSLQ: 1536 g33 + 8 Xs - 450 z7 + 180 z2 z5 + 45 z2^2 z3 = 0
g33_x = (-8*Xs + 450*z(7) - 180*z(2)*z(5) - 45*z(2)**2*z(3))/1536
print("gamma33 residual:", mp.nstr(g33_x-g33, 3))
# c7 = 40 g33 + (3/4) z3 z4 + 2339/512 z7  with z2^2 = (5/2) z4
a7 = Fr(375, 32)+Fr(2339, 512); a25 = -Fr(5*180, 192); a34 = -Fr(5*45, 192)*Fr(5, 2)+Fr(3, 4)
print("c7 = -(5/24) Xs + %s z7 + %s z2 z5 + %s z3 z4" % (a7, a25, a34))
c7 = -mp.mpf(5)/24*Xs + mp.mpf(a7.numerator)/a7.denominator*z(7) + mp.mpf(a25.numerator)/a25.denominator*z(2)*z(5) + mp.mpf(a34.numerator)/a34.denominator*z(3)*z(4)
c7g = 40*g33 + mp.mpf(3)/4*z(3)*z(4) + mp.mpf(2339)/512*z(7)
print("c7 closed form  :", mp.nstr(c7, 60), " diff to gamma formula", mp.nstr(c7-c7g, 3))
c7R = mp.mpf(open("sm6_c7_numeric.txt").read().strip())
print("c7 Richardson   :", mp.nstr(c7R, 50), " diff", mp.nstr(c7R-c7, 3))
# c6
Z51 = mp.quad(lambda x: (-mp.log(x))**4*mp.log(1-x)/(1+x), [0, mp.mpf(1)/2, 1], maxdegree=8)/24
c6 = (-173*z(6)+44*z(3)**2+155*z(5)*mp.log(2)-80*Z51)/32
c6R = mp.mpf("0.0429508278802840339294338048351352655102967276920578653499916")
print("c6 closed form  :", mp.nstr(c6, 60), " diff to Richardson", mp.nstr(c6-c6R, 3))
# PSLQ negative / positive tests for c7
mp.mp.dps = 47
print("c7 vs {z7, z2 z5, z3 z4}:", mp.pslq([c7R, z(7), z(2)*z(5), z(3)*z(4)], maxcoeff=10**8, maxsteps=10**6))
print("c7 vs {z7, z2 z5, z3 z4, Xs}:", mp.pslq([c7R, z(7), z(2)*z(5), z(3)*z(4), Xs], maxcoeff=10**5, maxsteps=10**6))
mp.mp.dps = 300
print("gamma33 vs {z7, z2 z5, z3 z4}:", mp.pslq([g33, z(7), z(2)*z(5), z(3)*z(4)], maxcoeff=10**12, maxsteps=10**6))
print("Xs vs {z7, z2 z5, z3 z4}:", mp.pslq([Xs, z(7), z(2)*z(5), z(3)*z(4)], maxcoeff=10**12, maxsteps=10**6))
print("Xs vs {z7, z2 z5, z3 z4, ln2*(wt 6 pure zeta), ...} skipped; see sm6_pslq7_output.txt for the 21-element basis")
# series check gamma32 = (1/16) sum_{n>=2} (-1)^n H_{n-1}/n^5
mp.mp.dps = 25
S = mp.nsum(lambda n: (-1)**int(n)*mp.harmonic(n-1)/n**5, [2, mp.inf])
g32 = mp.mpf("0.00164994679956980918331363168768428290424162265441038446989217")
print("gamma32 vs (1/16) sum:", mp.nstr(g32, 20), mp.nstr(S/16, 20))
