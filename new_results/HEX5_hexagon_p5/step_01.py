import sys, time
sys.path.insert(0, '../CX3k_review_draft8')
from polygon_engine import PolygonFamily, regular_vs
from mpmath import mp, mpf, log, pi, sqrt, zeta, nstr

t0 = time.time()
mp.dps = 70
vs, ws = regular_vs(3)   # n_i = (cos(i pi/3), sin(i pi/3)), i=1,2,3, weights 1
print("vs =", [tuple(nstr(x,8) for x in v) for v in vs])
PF = PolygonFamily(vs, ws, N=6)
print("wedges:", [(nstr(a,8), nstr(b,8)) for a,b in PF.wedges])
coeffs, maxdiff = PF.Pcoeffs(dps=60, kmax=3)
P3 = coeffs[3]
# paper closed form: P3(1) = (sqrt3/8)(16 log^3 2 + 3 pi^2 log2 - 96 log^2 2 + 96 log2 - 9 zeta(3) - 6 pi^2)
L2 = log(2)
P3_paper = sqrt(3)/8*(16*L2**3 + 3*pi**2*L2 - 96*L2**2 + 96*L2 - 9*zeta(3) - 6*pi**2)
print("Pd3(1) numeric =", nstr(P3, 55))
print("Pd3(1) paper   =", nstr(P3_paper, 55))
print("abs diff       =", nstr(abs(P3-P3_paper), 8))
print("rel diff       =", nstr(abs(P3-P3_paper)/abs(P3), 8))
print("engine maxdiff (last two TS levels, all k<=3) =", nstr(maxdiff, 8))
for k in range(4):
    print("P^(%d)(1) = %s" % (k, nstr(coeffs[k], 40)))
print("elapsed", time.time()-t0)
