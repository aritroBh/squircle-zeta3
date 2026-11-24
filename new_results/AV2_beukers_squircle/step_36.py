
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, nsum, inf
mp.dps = 110
exec(open('av2core.py').read())
z2, z3, ln2, ln3 = zeta(2), zeta(3), log(2), log(3)
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag

def sum_parity(p, want_odd):
    eps = mpf(1)/p
    S = mpf(0); n = 1
    tol = mpf(10)**(-(mp.dps-10))
    while True:
        if (n % 2 == 1) == want_odd:
            T = B2_term(n, eps)
            S += T
        else:
            T = B2_term(n, eps)
        if n > 40*p and fabs(T) < tol:
            break
        n += 1
    return S

S_odd = sum_parity(2, True)
S_even = sum_parity(2, False)
B2_2, _ = B2_series(2)
print("S_odd =", nstr(S_odd, 60))
print("S_even =", nstr(S_even, 60))
print("sum check:", nstr(S_odd+S_even-B2_2, 10))

def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**10):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if rel is None or rel[0] == 0:
        print("%s: (none)" % name); return
    resid = fabs(sum(mpf(rel[i])*([val]+basis)[i] for i in range(len(rel))))
    scale = max(fabs(mpf(rel[i])*([val]+basis)[i]) for i in range(len(rel)))
    print("%s: %s  resid=%s" % (name, rel, nstr(resid/scale, 5)))

cb = [mpf(1), z3, Cl2, pi*Cl2, pi, pi**2, ln3, ln2]
hunt("S_odd", S_odd, cb)
hunt("S_even", S_even, cb)
hunt("S_odd b", S_odd, [mpf(1), z3, pi, pi**2, ln2, ln2**2, pi*ln2])
hunt("S_even b", S_even, [mpf(1), z3, pi, pi**2, ln2, ln2**2, pi*ln2])
