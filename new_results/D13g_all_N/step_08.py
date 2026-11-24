
import math, time, importlib, sys
sys.path.insert(0, '.')
import step_37
importlib.reload(step_37)
from step_37 import certified_y

ys = lambda N: math.sqrt(math.pi/(2*N))
t0 = time.time()
r = certified_y(ys(800), ys(400), nsub=4)
print(f"N-box [400,800]: C = {r['C'].upper()}   ({time.time()-t0:.2f}s)")
print("   pieces*N^1.5:", {k: round(float(v.mid()),3) for k,v in r.items() if k!='C'})
