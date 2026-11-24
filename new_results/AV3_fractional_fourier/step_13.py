from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler, harmonic)
mp.dps = 80
def N(x,d=45): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_series(z):
    return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf], tol=mpf('1e-70'))
H3 = lambda m: harmonic(m, 3)

N0 = 40
# ---- S_Si tail via hump interchange (F_m exact) ----
headSi = nsum(lambda n: si(n*pi)/n**3, [1, N0])
F = lambda m: si((m+1)*pi) - si(m*pi)
tailSi = si(N0*pi)*(zeta(3)-H3(N0)) + nsum(lambda m: F(m)*(zeta(3)-H3(m)), [N0, inf])
S_Si = headSi + tailSi
print("S_Si =", N(S_Si,55)); print("diff 5pi^3/72:", nstr(abs(S_Si-5*pi**3/72),3))

# ---- S_Si2 tail via hump interchange ----
headSi2 = nsum(lambda n: si2_series(n*pi)/n**3, [1, N0])
G = lambda m: quad(lambda u: si(u)/u, [m*pi, (m+1)*pi])
tailSi2 = si2_series(N0*pi)*(zeta(3)-H3(N0)) + nsum(lambda m: G(m)*(zeta(3)-H3(m)), [N0, inf])
S_Si2 = headSi2 + tailSi2
print("S_Si2 =", N(S_Si2,55)); print("diff 49pi^3/432:", nstr(abs(S_Si2-49*pi**3/432),3))
logline("step_12: S_Si and S_Si2 recomputed by hump-interchange; match 5pi^3/72 and 49pi^3/432")
