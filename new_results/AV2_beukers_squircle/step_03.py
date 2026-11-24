
import json
from mpmath import mp, mpf, zeta, fabs, log
mp.dps = 110
exec(open('av2core.py').read())

print("=== validate B2_quad2 & B2_quad3 against series (110 dps) ===")
for p in [1, 2, 4, 16, 32, 64]:
    s, nn = B2_series(p)
    q2 = B2_quad2(p); q3 = B2_quad3(p)
    print("p=%3d  |s-q2|=%s  |s-q3|=%s" % (p, mp.nstr(fabs(s-q2),3), mp.nstr(fabs(s-q3),3)))

print("=== large p ===")
for p in [256, 4096, 65536, 2**20]:
    q = B2_quad2(p); q3 = B2_quad3(p)
    print("p=%7d B2=%s" % (p, mp.nstr(q, 42)))
    print("          B2-zeta3=%s   |q2-q3|=%s" % (mp.nstr(q-zeta(3),6), mp.nstr(fabs(q-q3),3)))

print("=== C2 direct ===")
mp.dps = 100
C2 = C2_direct()
print("C2 =", mp.nstr(C2, 55))
# cross-check via B2 at two large p: B2 - zeta3 ~ C2 eps^2 + zeta3/2 eps^3 + ...
mp.dps = 110
z3 = zeta(3); c4 = mpf(9)/4*zeta(4)
for p in [2**16, 2**20]:
    q = B2_quad2(p); eps = mpf(1)/p
    r = q - z3 - z3/2*eps**3 + c4*0 + mpf(3)/16*zeta(4)*eps**4
    print("p=2^%d: (B2-zeta3-zeta3/2 eps^3+3/16 z4 eps^4)/eps^2 = %s" % (int(log(p,2)), mp.nstr(r/eps**2, 20)))
