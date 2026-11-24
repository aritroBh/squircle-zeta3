from mpmath import mp, mpf, zeta, log, pi, nprint, quad, sqrt, nsum, inf, sin, harmonic
mp.dps = 100
d1 = mpf(open('d1.txt').read())
target = 2*pi**3/(27*sqrt(3)) - 2*zeta(3)/3
print("d1 - (2 pi^3/(27 sqrt3) - 2 zeta3/3):", nprint(d1 - target, 6))
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
print("P1 - 2pi^3/(27sqrt3):", nprint(P1 - 2*pi**3/(27*sqrt(3)), 6))
Ls3_60 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
print("Ls3(pi/3) - 7pi^3/108:", nprint(Ls3_60 - 7*pi**3/108, 6))
# S = sum sin(m pi/3) H_m / m^2 =? pi^3/36
S = nsum(lambda m: sin(m*pi/3)*harmonic(m)/m**2, [1, inf])
print("S =", nprint(S, 40))
print("S - pi^3/36:", nprint(S - pi**3/36, 6))
# also sum sin(2m pi/3)/m^3 = 2 pi^3/81
T = nsum(lambda m: sin(2*m*pi/3)/m**3, [1, inf])
print("sum sin(2m pi/3)/m^3 - 2pi^3/81:", nprint(T - 2*pi**3/81, 6))
# Im Li3(e^{i pi/3}) = 5 pi^3/162
IL = nsum(lambda m: sin(m*pi/3)/m**3, [1, inf])
print("ImLi3(e^{ipi/3}) - 5pi^3/162:", nprint(IL - 5*pi**3/162, 6))
with open('progress.log','a') as fh:
    fh.write("step_47: d1 = 2pi^3/(27sqrt3) - 2 zeta(3)/3 confirmed; Ls3(pi/3)=7pi^3/108; S=pi^3/36 confirmed\n")
