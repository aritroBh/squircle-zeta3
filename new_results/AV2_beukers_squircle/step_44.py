
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, quad, asin
mp.dps = 110
z2, z3, ln2, ln3 = zeta(2), zeta(3), log(2), log(3)
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag
B3_1 = mpf('0.24661730238492374056516806883959205599718610151827')
# P1 = int_0^1 -log s * 4 asin(s/2)/(s sqrt(4-s^2)) ds
P1 = quad(lambda s: -log(s)*4*asin(s/2)/(s*sqrt(4-s**2)) if s > 0 else 0, [0, 1])
print("P1 =", nstr(P1, 60))
P2 = B3_1 - P1
print("P2 = B3(1)-P1 =", nstr(P2, 60))
def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**9):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if rel is None or rel[0] == 0:
        print("%s: (none)" % name); return
    resid = fabs(sum(mpf(rel[i])*([val]+basis)[i] for i in range(len(rel))))
    scale = max(fabs(mpf(rel[i])*([val]+basis)[i]) for i in range(len(rel)))
    tag = "OK" if resid/scale < mpf('1e-40') else "WEAK"
    print("%s [%s]: %s  resid=%s" % (name, tag, rel, nstr(resid/scale, 5)))
hunt("P1", P1, [mpf(1), z2, z3, pi, pi**2, pi**3, ln2, ln2**2, ln2**3, Cl2, pi*Cl2, sqrt(3)*Cl2])
hunt("P1b", P1, [mpf(1), z3, pi*Cl2, pi**2, ln3, pi*ln3])
hunt("P2", P2, [mpf(1), z2, z3, pi, pi**2, pi**3, ln2, ln2**2, ln2**3, Cl2, pi*Cl2, sqrt(3)*Cl2])
hunt("P2b", P2, [mpf(1), z3, ln2, z2, pi, ln2**2, pi*ln2])
