
import math, time, sys
sys.path.insert(0,'.')
import step_37
from step_37 import certified_y
ys = lambda N: math.sqrt(math.pi/(2*N))

def try_cell(n1, n2, nsub=2):
    t0=time.time()
    try:
        r = certified_y(ys(n2), ys(n1), nsub=nsub)
        return float(r['C'].upper()), time.time()-t0
    except AssertionError as e:
        return f"FAIL: {e}", time.time()-t0

for (n1,n2) in [(400,600),(600,800),(400,500),(500,640),(640,800),(800,1200),(1200,1600)]:
    C,dt = try_cell(n1,n2)
    print(f"[{n1},{n2}] ratio {n2/n1:.3f}: C={C}  ({dt:.2f}s)")
