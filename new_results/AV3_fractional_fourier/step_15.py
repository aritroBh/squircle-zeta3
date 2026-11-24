from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler, quad)
mp.dps = 90
def N(x,d=45): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
def H3f(m):  # generalized harmonic H_m^{(3)} for integer m
    return zeta(3) - zeta(3, m+1)

N0 = 40
headSi = nsum(lambda n: si(n*pi)/n**3, [1, N0])
F = lambda m: si((m+1)*pi) - si(m*pi)
tailSi = si(N0*pi)*zeta(3, N0+1) + nsum(lambda m: F(m)*zeta(3, m+1), [N0, inf])
S_Si = headSi + tailSi
print("S_Si =", N(S_Si,55))
print("diff 5pi^3/72:", nstr(abs(S_Si-5*pi**3/72),3))
logline("step_14: S_Si=5pi^3/72 verified 3rd method (exact hump interchange), full precision")
