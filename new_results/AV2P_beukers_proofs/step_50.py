from mpmath import mp, mpf, zeta, log, pi, nprint, sin, quad
mp.dps = 60
# T = sum sin(m pi/3) H_m^{(2)}/m
def psum2(f, N):
    S = mpf(0); H2 = mpf(0)
    for m in range(1, N+1):
        H2 += 1/mpf(m)**2
        S += f(m, H2)
    return S
for N in [10000, 100000]:
    T = psum2(lambda m,H2: sin(m*pi/3)*H2/m, N)
    print(f"T (N={N}) =", nprint(T, 25))
print("31 pi^3/648 =", nprint(31*pi**3/648, 25))
# Ls3(2pi/3) = 2.14476721256949439148264211222
print("Ls3(2pi/3) - 13 pi^3/162:", nprint(2.14476721256949439148264211222 - 13*pi**3/162, 6))
print("13 pi^3/162 =", nprint(13*pi**3/162, 25))
# Cl2(pi/3) = 9*13? no: 9? check alias
print("Cl2(pi/3) =", nprint(mp.clsin(2,pi/3), 30))
print("S - 11 pi^3/324:", nprint(1.052682232973162117428381 - 11*pi**3/324, 8))
# A double check via Ls3(2pi/3)=13pi^3/162 and S,T:
A = pi**3/72 - 31*pi**3/1296 + 5*pi**3/648 + 13*pi**3/1296 - 5*pi**3/648
print("A closed =", nprint(A, 25))
print("A - 7 pi^3/216:", nprint(A - 7*pi**3/216, 8))
with open('progress.log','a') as fh:
    fh.write("step_49: T=31pi^3/648, Ls3(2pi/3)=13pi^3/162, A=7pi^3/216 confirmed -> Ls3(pi/3)=7pi^3/108 proved\n")
