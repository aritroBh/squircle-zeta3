step08 = r'''import time, importlib.util
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
'''
open('step_08.py','w').write(step08)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_08.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-300:])
