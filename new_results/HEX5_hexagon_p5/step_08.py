import time, importlib.util
spec = importlib.util.spec_from_file_location('polygon_engine', '../CX3k_review_draft8/polygon_engine.py')
PE = importlib.util.module_from_spec(spec); spec.loader.exec_module(PE)
from mpmath import mp, nstr
t0 = time.time()
vs, ws = PE.regular_vs(3)
PF = PE.PolygonFamily(vs, ws, N=7)
coeffs, maxdiff = PF.Pcoeffs(dps=50, kmax=5)
print("original engine, dps=50, kmax=5, elapsed %.1fs, maxdiff=%s" % (time.time()-t0, nstr(maxdiff,6)))
for k in range(6):
    print("P^(%d)(1) = %s" % (k, nstr(coeffs[k], 45)))
