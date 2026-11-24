from mpmath import mp, mpf, zeta, log, pi, nprint, quad, binomial, nsum, inf, sin, harmonic, sqrt
mp.dps = 60
Ls3 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
A = quad(lambda th: log(2*mp.sin(th))**2, [0, pi/6])
print("Ls3(pi/3) =", nprint(Ls3, 30))
print("2*A =", nprint(2*A, 30))
S4 = mpf(0)
for m in range(0, 400):
    S4 += binomial(2*m,m)/(mpf(16)**m*(2*m+1)**3)
print("S4 = sum C/(16^m(2m+1)^3) =", nprint(S4, 30))
print("Ls3/2 =", nprint(Ls3/2, 30))
print("Ls3 - 2 S4:", nprint(Ls3 - 2*S4, 5))
# so Ls3 = 2*sumC/(16^m(2m+1)^3): hence sum = Ls3/2 = 7 pi^3/216
print("S4 - 7pi^3/216:", nprint(S4 - 7*pi**3/216, 5))
# direct central sum S5 = sum C(1/2)^{2m+1}/(2m+1)^3
S5 = mpf(0)
for m in range(0, 400):
    S5 += binomial(2*m,m)*mpf('0.5')**(2*m+1)/(2*m+1)**3
print("S5 =", nprint(S5, 30))
# sin harmonic sums, direct partial sums
S = mpf(0); Hk = mpf(0)
for m in range(1, 20001):
    Hk += 1/mpf(m)
    S += sin(m*pi/3)*Hk/m**2
print("S (direct, m<=20000) =", nprint(S, 30))
print("S - pi^3/36:", nprint(S - pi**3/36, 5))
