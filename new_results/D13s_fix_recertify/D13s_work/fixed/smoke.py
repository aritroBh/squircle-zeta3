
import math, time
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
for (n1,n2) in [(100.0,105.0),(390.0,400.0),(400.0,600.0),(1e6,1.5e6)]:
    t0=time.time()
    r = certified_y38(ys(n2), ys(n1), J=24, nsub=4)
    print(f"[{n1},{n2}] C = {float(r[chr(67)].upper()):.4f}  ({time.time()-t0:.1f}s)")
