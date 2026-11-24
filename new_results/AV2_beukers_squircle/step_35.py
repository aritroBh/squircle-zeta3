
from mpmath import mp, mpf, zeta, fabs, log, pi, sqrt, nstr, pslq, polylog, nsum, inf
mp.dps = 105
from math import comb
z2, z3, ln2, ln3 = zeta(2), zeta(3), log(2), log(3)
Cl2 = (polylog(2, mp.e**(1j*pi/3))).imag
def H(m):
    return nsum(lambda k: 1/k, [1, int(m)])
SA = nsum(lambda k: (H(2*k)-H(k))/(int(k)**2*comb(2*int(k), int(k))), [1, inf])
SB = nsum(lambda k: 1/(int(k)**3*comb(2*int(k), int(k))), [1, inf])
SC = nsum(lambda k: 1/(int(k)**4*comb(2*int(k), int(k))), [1, inf])
print("SA =", nstr(SA, 55))
print("SB =", nstr(SB, 55))
print("SC =", nstr(SC, 55))
def hunt(name, val, basis, tol=mpf('1e-45'), maxcoeff=10**10):
    rel = pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    ok = rel is not None and rel[0] != 0
    print("%s: %s" % (name, rel if ok else "(none %s)" % rel))
cb = [mpf(1), z3, Cl2, pi*Cl2, pi, pi**2, ln3, sqrt(3), ln2]
hunt("SA", SA, cb)
hunt("SB", SB, cb)
hunt("SC", SC, cb)
hunt("SA2", SA, [mpf(1), z3, pi*Cl2, pi, ln3, ln2, pi**2])
hunt("SB2", SB, [mpf(1), z3, pi*Cl2, pi, ln3, ln2, pi**2])
# B2(2) even part check
exec(open('av2core.py').read())
B2_2, _ = B2_series(2)
S1  = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * H(m)/(2*m+1), [1, inf])
S1p = nsum(lambda m: mpf(comb(2*int(m), int(m)))/16**int(m) * H(2*m)/(2*m+1), [1, inf])
S2s = Cl2 - 1
S3s = (pi-3)/3
odd = pi/8*((S1-S1p) + S2s + 2*ln2*S3s)
print("odd part =", nstr(odd, 40))
print("odd+ (SA+SB)/8 =", nstr(odd + (SA+SB)/8, 40))
print("B2(2) =", nstr(B2_2, 40))
print("diff:", nstr(B2_2 - (odd+(SA+SB)/8), 10))
