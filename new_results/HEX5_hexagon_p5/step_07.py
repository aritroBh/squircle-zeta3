import re
from mpmath import mp, mpf, log, pi, zeta, sqrt, nstr, pslq, fabs

txt = open('p5_run150.txt').read()
val_str = re.search(r'P5/sqrt3 = ([0-9.]+)', txt).group(1)
REL = [128, 0, -30720, 4920, 61440, 0, -6240, -30720, -448, -5760, 1920, 5120, -11340, -5580, 123, 1440, -160, -256]
names = ['1','log2','pi^2','log2^2','z3','pi^2 log2','log2^3','pi^4','z3 log2',
         'pi^2 log2^2','log2^4','z5','z2 z3','pi^4 log2','z3 log2^2','pi^2 log2^3','log2^5']

mp.dps = 170
L2 = log(2)
basis = [mpf(1), L2, pi**2, L2**2, zeta(3), pi**2*L2, L2**3, pi**4, zeta(3)*L2,
         pi**2*L2**2, L2**4, zeta(5), (pi**2/6)*zeta(3), pi**4*L2, zeta(3)*L2**2,
         pi**2*L2**3, L2**5]
val = mpf(val_str)
x = [val] + basis
resid = fabs(sum(mpf(c)*xx for c, xx in zip(REL, x)))
print("independent check at dps=170: residual of 128-relation =", nstr(resid, 6))
print("max|coeff| =", max(abs(c) for c in REL))

# reconstruct value from relation, compare digits
recon = -sum(mpf(REL[i+1])*basis[i] for i in range(17))/mpf(REL[0])
diff = fabs(recon - val)
print("reconstruction |recon - stored| =", nstr(diff, 6))
import mpmath
nd = -mpmath.log10(diff/fabs(val))
print("relative agreement digits ~", nstr(nd, 6))

# zeta coefficients
print("z5 coeff    = 11340/128 =", nstr(mpf(11340)/128, 12), " ; claimed 2835/32 =", nstr(mpf(2835)/32, 12), "; equal:", mpf(11340)/128 == mpf(2835)/32)
print("z2z3 coeff  = 5580/128  =", nstr(mpf(5580)/128, 12), " ; claimed 1395/32 =", nstr(mpf(1395)/32, 12), "; equal:", mpf(5580)/128 == mpf(1395)/32)

# subtraction variant: fix the two zeta terms, PSLQ the rest (dim 16)
for DPS in [100, 125]:
    mp.dps = DPS
    L2 = log(2)
    b5 = zeta(5); b23 = (pi**2/6)*zeta(3)
    R = mpf(val_str) - mpf(2835)/32*b5 - mpf(1395)/32*b23
    keep = [0,1,2,3,4,5,6,7,8,9,10,13,14,15,16]
    L2v = log(2)
    bb = [mpf(1), L2v, pi**2, L2v**2, zeta(3), pi**2*L2v, L2v**3, pi**4, zeta(3)*L2v,
          pi**2*L2v**2, L2v**4, pi**4*L2v, zeta(3)*L2v**2, pi**2*L2v**3, L2v**5]
    xx = [R] + bb
    rel = pslq(xx, tol=mpf(10)**(-(DPS-25)), maxcoeff=10**9, maxsteps=10000)
    resid2 = fabs(sum(mpf(rel[i])*xx[i] for i in range(len(xx)))) if rel else None
    print("subtraction PSLQ DPS=%d: residual=%s" % (DPS, nstr(resid2,6) if resid2 else 'None'))
    print("  rel =", rel)
    exp = [128, 0, -30720, 4920, 61440, 0, -6240, -30720, -448, -5760, 1920, 5120, 123, 1440, -160, -256]
    print("  matches 128-relation restricted:", rel == exp)
