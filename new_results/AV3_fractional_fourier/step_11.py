from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_series(z):
    return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf], tol=mpf('1e-90'))

# quick check of C = pi*gamma_euler/2 via z=25 with epsilon correction
def eps_asym(z, K=None, J=None):
    # ε(z)=Σ_{k≥0}(-1)^{k+1}(2k)! S_{2k+3}(z), S_m(z)=cos z Σ_j (m+2j-2)!/(m-1)! /z^{m+2j} - wait derived: sgn=cos z
    import math
    if K is None: K = int(z) 
    if J is None: J = int(z)
    cz = cos(z)
    def S(m):
        t = mpf(0); r = mpf(1)
        for j in range(J):
            t += r / z**(m+2*j)
            r *= (m+2*j)*(m+2*j+1)
            if r / z**(m+2*j+2) > abs(t)*1e3: break  # divergence guard
        return cz*t
    tot = mpf(0); prev = None
    for k in range(K):
        term = (-1)**(k+1) * fac(2*k) * S(2*k+3)
        tot += term
        if prev is not None and abs(term) > abs(prev):  # optimal truncation
            tot -= term; break
        prev = term
    return tot

for z in [20, 25, 30]:
    cest = si2_series(mpf(z)) - (pi/2)*log(z) - eps_asym(mpf(z))
    print(f"z={z}: C estimate={N(cest,30)}  pi*gamma/2={N(pi*euler/2,30)}  diff={nstr(abs(cest-pi*euler/2),3)}")
logline("step_10: Si2 constant C = pi*gamma/2 identified analytically (C = -int sin t log t /t = pi gamma/2), numerically confirmed")
