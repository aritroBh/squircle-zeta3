from mpmath import mp, mpf, log, pi, nprint, quad, clsin
mp.dps = 60
# check int_0^{pi/6} logcos t dt vs -(pi/6)log2 + Cl2(2pi/3)/2
I1 = quad(lambda t: log(mp.cos(t)), [0, pi/6])
c1 = -(pi/6)*log(2) + clsin(2, 2*pi/3)/2
print("I1 =", nprint(I1, 30), "\nc1 =", nprint(c1, 30), "\ndiff:", nprint(I1-c1, 5))
# check int_0^{pi/12} logcos u du vs -(pi/12)log2 + Cl2(5pi/6)/2
I2 = quad(lambda t: log(mp.cos(t)), [0, pi/12])
c2 = -(pi/12)*log(2) + clsin(2, 5*pi/6)/2
print("I2 =", nprint(I2, 30), "\nc2 =", nprint(c2, 30), "\ndiff:", nprint(I2-c2, 5))
# check the integral formula int_0^th logcos x dx = -th log2 + Cl2(pi-2th)/2
th = mpf('0.4')
lhs = quad(lambda x: log(mp.cos(x)), [0, th])
rhs = -th*log(2) + clsin(2, pi-2*th)/2
print("formula diff:", nprint(lhs-rhs, 5))
# Cl2(5pi/6) check vs identity
print("Cl2(5pi/6) =", nprint(clsin(2,5*pi/6), 30))
print("Cl2(pi/6)-Cl2(pi/3)/2 =", nprint(clsin(2,pi/6)-clsin(2,pi/3)/2, 30))
