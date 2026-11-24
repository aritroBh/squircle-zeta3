
import mpmath as mp
mp.mp.dps = 90

z3 = mp.zeta(3); ln2 = mp.log(2); pi2 = mp.pi**2; sq3 = mp.sqrt(3)
P3 = mp.mpf("-5.14545501713561583717168113112397507111555356641974172092761")
A3 = P3/6
c1 = 2*sq3*ln2
Ig1 = pi2*sq3/48
c2 = ln2**2*sq3 + 3*Ig1
Ig1sq_half = sq3*z3/64
Ialpha_s_g1 = 3*sq3*z3/64
Ig2 = (A3 - c1 + 2*c2 - (ln2**3/3)*sq3)/3 - ln2*Ig1 - Ig1sq_half - Ialpha_s_g1
print("backed-out int g2 =", mp.nstr(Ig2, 60))
print("-sq3 z3/8 =", mp.nstr(-sq3*z3/8, 60))
basis = [z3, z3*ln2, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
M = [Ig2/sq3] + basis
rel = mp.pslq(M, tol=mp.mpf('1e-80'), maxcoeff=10**12, maxsteps=10000)
print("PSLQ Ig2/sq3:", rel)
if rel:
    print(" resid:", mp.nstr(mp.fsum(mp.mpf(r)*x for r,x in zip(rel,M)),8))
    print(" zeta(3) coefficient:", mp.mpf(-rel[1])/rel[0])
