
import importlib, sys, math, time
import step_38; importlib.reload(step_38)
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
r = certified_y38(ys(150), ys(100), nsub=4)
print("[100,150] pieces*N^1.5:", {k: round(float(v.mid()),2) for k,v in r.items() if k!='C'}, "C=", float(r['C'].upper()))
# probe dyadic cells at large N
for (n1,n2) in [(400,600),(600,900),(10**4,1.5*10**4),(10**6,1.5*10**6),(10**8,1.5*10**8)]:
    n2=int(n2)
    t0=time.time()
    try:
        r = certified_y38(ys(n2), ys(n1), nsub=2)
        print(f"[{n1},{n2}] r={n2/n1:.2f}: C={float(r['C'].upper()):.1f}  ({time.time()-t0:.2f}s) TLf={float(r['TLf'].mid()):.1f} W3={float(r['W3'].mid()):.1f} lamok")
    except AssertionError as e:
        print(f"[{n1},{n2}]: FAIL {e}  ({time.time()-t0:.2f}s)")
