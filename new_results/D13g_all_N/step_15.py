
src = open('step_38.py').read()
src = src.replace("r = certified_y(", "r = certified_y38(")
open('step_38.py','w').write(src)
import importlib, sys, math, time
sys.path.insert(0,'.')
import step_38; importlib.reload(step_38)
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
# regression: the four D13f cells
for (n1,n2) in [(100,150),(150,225),(225,340),(340,400)]:
    t0=time.time()
    r = certified_y38(ys(n2), ys(n1), nsub=4)
    print(f"[{n1},{n2}]: C={float(r['C'].upper()):.2f}  TLf*N1.5={float(r['TLf'].mid()):.3f}  ({time.time()-t0:.2f}s)")
