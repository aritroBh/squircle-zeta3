from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler)
import json, os
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

hc = {int(r): mpf(v.replace('pi','').strip()) for r,v in {}.items()}
raw = json.load(open("h_coeffs.json"))
import re
def ev(s):
    return mpf(s.replace('pi', f'*{mp.nstr(pi,80)}').replace('**','^')) if False else None
# safer: rebuild h with mpmath from sympy-free recursion
mp.dps = 100
pmax = 30
CA = [mpf(0)]*pmax; CB = [mpf(0)]*pmax
CA[0] = mpf(0); CB[0] = mpf(2)
for j in range(1, pmax):
    CA[j] = -j*CB[j-1]
    CB[j] = pi**j + j*CA[j-1]
KMAX = 12
h = {}
for k in range(KMAX):
    for (P, moments, fk) in [(2*k+2, CA, fac(2*k)), (2*k+3, CB, fac(2*k+1))]:
        for j in range(0, pmax - P):
            r = P + j
            if r > pmax: continue
            coeff = (-1)**j * gamma(P+j)/gamma(P)/fac(j) * moments[j]
            h[r] = h.get(r, mpf(0)) + (-1)*(-1)**k * fk * coeff
print("h check:", nstr(h[3],6), nstr(h[4],6), nstr(h[5],6))  # expect -6, 9pi, 100-18pi^2

def H_asym(m):  # (-1)^m sum h_r (m pi)^{-r}, r=3..26
    s = mpf(0); z = m*pi
    for r in range(3, 27):
        s += h[r]/z**r
    return (-1)**int(m % 2) * s

# exact humps m=1..200 at elevated precision
if os.path.exists("humps.json"):
    humps = {int(k): mpf(v) for k,v in json.load(open("humps.json")).items()}
else:
    mp.dps = 260
    def si2_series(z):
        return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf], tol=mpf('1e-240'))
    humps = {}
    for m in range(1, 201):
        humps[m] = si2_series((m+1)*pi) - si2_series(m*pi) - (pi/2)*log((m+1)/m)
    json.dump({str(k): mp.nstr(v, 240) for k,v in humps.items()}, open("humps.json","w"))
    mp.dps = 100
print("humps ready:", len(humps))
json.load(open("progress.log")) if False else None
logline("step_16: h_r coefficients rebuilt; exact humps m<=200 cached")
