import sys
sys.path.insert(0, '../AV2_beukers_squircle')
from mpmath import mp, mpf, nsum, inf, pi, zeta, log, sin, cos, sqrt, clsin, clcos, nprint

mp.dps = 50

# Classical constants, computed independently
Cl2_pi3 = clsin(2, pi/3)      # Cl_2(pi/3)
Cl2_pi6 = clsin(2, pi/6)
Cl3_pi3 = clcos(3, pi/3)      # Cl_3(pi/3) = sum cos(k pi/3)/k^3
z3 = zeta(3)

print("Cl2(pi/3) =", nprint(Cl2_pi3, 40))
print("Cl2(pi/6) =", nprint(Cl2_pi6, 40))
print("Cl3(pi/3) =", nprint(Cl3_pi3, 40))
print("zeta(3)/3 =", nprint(z3/3, 40))
print("Cl3(pi/3) - zeta(3)/3 =", nprint(Cl3_pi3 - z3/3, 10))

# U = sum 1/(n^3 C(2n,n))
from mpmath import binomial
U = nsum(lambda n: 1/(mpf(n)**3 * binomial(2*n, n)), [1, inf])
print("U =", nprint(U, 40))
cand = (2*pi/3)*Cl2_pi3 - 4*z3/3
print("(2pi/3)Cl2(pi/3) - 4 zeta(3)/3 =", nprint(cand, 40))
print("U - cand =", nprint(U - cand, 10))
cand2 = -(2*pi/3)*Cl2_pi3 + 2*z3  # another guess
print("U - cand2 =", nprint(U - cand2, 10))

# T = sum (H_{2n}-H_n)/(n^2 C(2n,n))
from mpmath import harmonic
T = nsum(lambda n: (harmonic(2*n)-harmonic(n))/(mpf(n)**2*binomial(2*n,n)), [1, inf])
print("T =", nprint(T, 40))
# S = B2(1) = sum (H_{2n}-H_{n-1})/(n^2 C(2n,n))
S = nsum(lambda n: (harmonic(2*n)-harmonic(n-1))/(mpf(n)**2*binomial(2*n,n)), [1, inf])
print("S = B2(1) series =", nprint(S, 40))
print("2 zeta(3)/3 =", nprint(2*z3/3, 40))
print("S - 2zeta3/3 =", nprint(S - 2*z3/3, 10))
print("T + U =", nprint(T+U, 40))
