
from mpmath import mp, mpf, zeta, fabs, log, log1p, pi, sqrt, sin, nstr, pslq, polylog, euler
mp.dps = 110
exec(open('av2core.py').read())

# values (>=100 digits)
B2_1, _ = B2_series(1)           # 0.801371268773062856933...
B3_1 = mpf('0.24661730238492374056516806883959205599718610151827')
B2_2, _ = B2_series(2)           # 1.059368335799353241...
B2_4, _ = B2_series(4)
print("B2(1) =", nstr(B2_1, 60))
print("B3(1) =", nstr(B3_1, 60))
print("B2(2) =", nstr(B2_2, 60))

# constants
z2, z3, z4, z5 = zeta(2), zeta(3), zeta(4), zeta(5)
ln2, ln3, gam = log(2), log(3), euler
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag   # Clausen Cl_2(pi/3)
print("Cl2(pi/3) =", nstr(Cl2, 30))
pis3 = pi/sqrt(3)

def hunt(name, val, basis, tol=mpf('1e-60'), maxcoeff=10**10):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    # filter trivial: require rel[0] != 0
    ok = rel is not None and rel[0] != 0
    print("%s: %s" % (name, rel if ok else ("(none; got %s)" % (rel,))))
    return rel

print("=== B2(1) hunts ===")
hunt("z3-pi2", B2_1, [mpf(1), z2, z3])
hunt("z3-claus", B2_1, [mpf(1), z2, z3, pi*Cl2, sqrt(3)*Cl2, pis3, pi*pis3])
hunt("z3-pi-ln", B2_1, [mpf(1), z2, z3, ln2, ln3, z2*ln2, z2*ln3, ln2**2, ln3**2, ln2*ln3])
print("=== B3(1) hunts ===")
hunt("z3-pi2", B3_1, [mpf(1), z2, z3])
hunt("z3-claus", B3_1, [mpf(1), z2, z3, pi*Cl2, sqrt(3)*Cl2, pis3, pi*pis3])
hunt("z3-pi-ln", B3_1, [mpf(1), z2, z3, ln2, ln3, z2*ln2, z2*ln3, ln2**2, ln3**2, ln2*ln3])
hunt("frac-z3", B3_1/z3, [mpf(1)])
print("=== B2(2): Q + Q ln2 + pi Q + pi ln2 Q ===")
hunt("q-ln2-pi", B2_2, [mpf(1), ln2, pi, pi*ln2], maxcoeff=10**12)
hunt("q-ln2-pi-ext", B2_2, [mpf(1), ln2, pi, pi*ln2, z2, z3, ln2**2, pi*ln2**2, z2*ln2])
print("=== B2(4) ===")
hunt("b2-4", B2_4, [mpf(1), ln2, pi, pi*ln2, z2, z3, ln2**2, pi*ln2**2, z2*ln2, sqrt(2), log(1+sqrt(2)), pi*sqrt(2)])
