from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, catalan, binomial, polylog
mp.dps = 80
G = catalan; Cl3 = clsin(2, pi/3); Cl23 = clsin(2, 2*pi/3)
# U = sum C(2m,m) H_m/(16^m m) =? 2 Cl_2(2pi/3)
K = 400
U = mpf(0); Hk = mpf(0)
for m in range(1, K):
    Hk += 1/mpf(m)
    U += binomial(2*m,m)*Hk/(mpf(16)**m*m)
print("U =", nprint(U, 40))
print("2 Cl2(2pi/3) =", nprint(2*Cl23, 40))
print("U - 2Cl2(2pi/3):", nprint(U-2*Cl23, 6))
print("Cl2(2pi/3)/Cl2(pi/3) =", nprint(Cl23/Cl3, 30))

# B2(2) pieces A, B
A = mpf('0.5')*quad(lambda a: polylog(2, mp.sin(a))/mp.sin(a), [0, pi/4])
B = mpf('0.5')*quad(lambda a: log(mp.sin(a))*log(1-mp.sin(a))/mp.sin(a), [0, pi/4])
print("\nA =", nprint(A, 40))
print("B =", nprint(B, 40))
print("A+B =", nprint(A+B, 40))
print("zeta3/12 + pi G/3 =", nprint(zeta(3)/12 + pi*G/3, 40))
print("A+B - target:", nprint(A+B - zeta(3)/12 - pi*G/3, 6))
