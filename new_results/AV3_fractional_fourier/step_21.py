from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler)
import json
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
mp.dps = 100
def N(x,d=42): return nstr(x,d)

# rebuild h_r
pmax = 30
CA = [mpf(0)]*pmax; CB = [mpf(0)]*pmax
CA[0] = mpf(0); CB[0] = mpf(2)
for j in range(1, pmax):
    CA[j] = -j*CB[j-1]; CB[j] = pi**j + j*CA[j-1]
KMAX = 12
h = {}
for k in range(KMAX):
    for (P, moments, fk) in [(2*k+2, CA, fac(2*k)), (2*k+3, CB, fac(2*k+1))]:
        for j in range(0, pmax - P):
            r = P + j
            if r > pmax: continue
            coeff = (-1)**j * gamma(P+j)/gamma(P)/fac(j) * moments[j]
            h[r] = h.get(r, mpf(0)) + (-1)*(-1)**k * fk * coeff

humps = {int(k): mpf(v) for k,v in json.load(open("humps.json")).items()}
def Hfun(m):  # valid for real m >= 1
    if m <= 120:
        # exact: interpolate? no — need exact at real m too; nsum only samples integers for alternating? not guaranteed
        return None
    s = mpf(0); z = m*pi
    for r in range(3, 27): s += h[r]/z**r
    return cos(m*pi) * s

# nsum for alternating series samples at integer points (it uses even/odd splitting), so integer m only needed for extrapolation? 
# Safer: sum E = -sum_{m>=1} H_m H3_m directly with manual splitting:
#   m=1..120 exact from cache; m=121..M asymptotic; tail via nsum on asymptotic from M+1.
def H3f(m): return zeta(3) - zeta(3, m+1)
E_head = -sum(humps[m]*H3f(m) for m in range(1,121))
M = 400
def H_asym_int(m):
    s = mpf(0); z = m*pi
    for r in range(3, 27): s += h[r]/z**r
    return ((-1)**m) * s
E_mid = -nsum(lambda m: H_asym_int(m)*H3f(m), [121, M])
E_tail = -nsum(lambda m: H_asym_int(int(m))*H3f(m), [M+1, inf])
E = E_head + E_mid + E_tail
print("E =", N(E,45))
S_Si2 = (pi/2)*(log(pi)+euler)*zeta(3) - (pi/2)*zeta(3, derivative=1) + E
print("S_Si2 =", N(S_Si2,45))
print("49pi^3/432 =", N(49*pi**3/432,45))
print("diff:", nstr(abs(S_Si2 - 49*pi**3/432), 3))
logline("step_19: S_Si2 computed via E decomposition; compared with 49pi^3/432")
