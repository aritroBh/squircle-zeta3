
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, euler, nsum, inf, quad, sin, cos
mp.dps = 105
exec(open('av2core.py').read())
z2, z3, ln2 = zeta(2), zeta(3), log(2)
from math import comb, factorial
from fractions import Fraction

# --- B2(2) decomposition ---
# odd part = 3/32 + zeta(3)/12 ; even part E2 = (pi/8)(S1 - S2 + 2 ln2 S3)
# S1 = sum_{m>=1} C(2m,m) H_m / (16^m (2m+1)), S2 = sum C/(16^m (2m+1)^2), S3 = sum C/(16^m (2m+1))
def H(n):
    return sum(Fraction(1, k) for k in range(1, n+1))
S1 = mpf(0); S2 = mpf(0); S3v = mpf(0)
for m in range(1, 400):
    c = mpf(comb(2*m, m)) / 16**m
    Hm = nsum(lambda k: 1/k, [1, m])
    S1 += c*Hm/(2*m+1)
    S2 += c/(2*m+1)**2
    S3v += c/(2*m+1)
# tail ~ 1/sqrt(m): add Euler-Maclaurin-ish tail estimate via nsum extrapolation
S1t = S1 + nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * (log(m)+euler)/(2*m+1), [400, inf]) if False else None
print("S1 (m<400) =", nstr(S1, 60))
print("S2 (m<400) =", nstr(S2, 60))
print("S3 (m<400) =", nstr(S3v, 60))
print("tails ~1e-4/sqrt... need acceleration - use nsum with Richardson instead:")
S1f = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * nsum(lambda k: 1/k, [1, int(m)])/(2*m+1), [1, inf])
S2f = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m)/(2*m+1)**2, [1, inf])
S3f = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m)/(2*m+1), [1, inf])
print("S1 =", nstr(S1f, 60))
print("S2 =", nstr(S2f, 60))
print("S3 =", nstr(S3f, 60))
B2_2, _ = B2_series(2)
E2 = B2_2 - (mpf(3)/32 + z3/12)
print("E2 = B2(2) - 3/32 - zeta3/12 =", nstr(E2, 60))
print("(pi/8)(S1-S2+2ln2 S3) =", nstr(pi/8*(S1f-S2f+2*ln2*S3f), 60))
print("diff:", nstr(E2 - pi/8*(S1f-S2f+2*ln2*S3f), 10))
print("=== hunts for S1, S2, S3, E2 ===")
def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**10):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    ok = rel is not None and rel[0] != 0
    print("%s: %s" % (name, rel if ok else "(none %s)" % rel))
hunt("S2", S2f, [mpf(1), pi, pi**2, pi**3])
hunt("S3", S3f, [mpf(1), pi, pi**2, pi**3])
hunt("S1", S1f, [mpf(1), pi, pi**2, ln2, ln2**2, pi*ln2, pi**3, pi*ln2**2])
hunt("E2", E2, [mpf(1), pi, pi**2, ln2, ln2**2, pi*ln2, pi**3, pi*ln2**2, z3])
