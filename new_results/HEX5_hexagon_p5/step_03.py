import sys, time
sys.path.insert(0, '../CX3k_review_draft8')
from polygon_engine import PolygonFamily, regular_vs
from mpmath import mp, mpf, pi, nstr

mp.dps = 40
vs, ws = regular_vs(3)
PF = PolygonFamily(vs, ws, N=7)
ths = [0.7, 1.1, 1.9, 2.4, 3.0, 3.4]
print("check pi/3 periodicity of integrand series:")
base = None
for th in ths:
    F = PF.integrand_series(mpf(th))
    Fp = PF.integrand_series(mpf(th) + pi/3)
    d = max(abs(a-b) for a,b in zip(F,Fp))
    print(" th=%s  max|F(th)-F(th+pi/3)| = %s" % (th, nstr(d,4)))
print("check symmetry about wedge midpoint pi/3 (=60 deg):")
mid = pi/3
for s in [0.1, 0.3, 0.45]:
    F1 = PF.integrand_series(mid + mpf(s))
    F2 = PF.integrand_series(mid - mpf(s))
    d = max(abs(a-b) for a,b in zip(F1,F2))
    print(" s=%s  max|F(mid+s)-F(mid-s)| = %s" % (s, nstr(d,4)))
# are the three wedge integrals equal? compare series at corresponding points
a1,b1 = PF.wedges[0]; a2,b2 = PF.wedges[1]
print("wedge widths:", nstr(b1-a1,10), nstr(b2-a2,10))
