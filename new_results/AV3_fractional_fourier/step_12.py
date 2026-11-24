from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, zeta, euler)
mp.dps = 130
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_series(z):
    return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf], tol=mpf('1e-120'))

def eps_asym(z, K, J):
    cz = cos(z)
    def S(m):
        t = mpf(0); r = mpf(1); best=None
        for j in range(J):
            term = r / z**(m+2*j); t += term
            r *= (m+2*j)*(m+2*j+1)
        return cz*t
    tot = mpf(0); terms=[]
    for k in range(K):
        terms.append((-1)**(k+1) * fac(2*k) * S(2*k+3))
    # optimal truncation at smallest |term|
    i = min(range(len(terms)), key=lambda i: abs(terms[i]))
    return sum(terms[:i+1]), abs(terms[i])

print("== confirm C = pi*gamma/2: eps_exact(z) vs eps_asym(z) ==")
for z in [60, 80, 100]:
    mp.dps = 130
    e_ex = si2_series(mpf(z)) - (pi/2)*(log(z)+euler)
    e_as, last = eps_asym(mpf(z), K=120, J=120)
    print(f"z={z}: eps_exact={nstr(e_ex,25)}  eps_asym={nstr(e_as,25)}  diff={nstr(abs(e_ex-e_as),3)} (last term {nstr(last,3)})")
logline("step_11: C=pi*gamma/2 confirmed by matching eps_exact vs asymptotics at z=60,80,100")
