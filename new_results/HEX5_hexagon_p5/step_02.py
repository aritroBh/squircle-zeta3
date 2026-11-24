import sys, time
sys.path.insert(0, '../CX3k_review_draft8')
from polygon_engine import PolygonFamily, regular_vs
from mpmath import mp, mpf, log, pi, sqrt, zeta, nstr

DPS = int(sys.argv[1]); TAG = sys.argv[2]
t0 = time.time()
mp.dps = DPS + 15
vs, ws = regular_vs(3)
PF = PolygonFamily(vs, ws, N=7)
coeffs, maxdiff = PF.Pcoeffs(dps=DPS, kmax=5)
out = []
out.append("DPS=%s  elapsed=%.1fs  engine maxdiff(last 2 TS levels, k<=5)=%s" % (DPS, time.time()-t0, nstr(maxdiff,6)))
for k in range(6):
    out.append("P^(%d)(1) = %s" % (k, nstr(coeffs[k], DPS)))
out.append("P5/sqrt3 = " + nstr(coeffs[5]/sqrt(3), DPS))
txt = "\n".join(out)
print(txt)
open(TAG,'w').write(txt+"\n")
