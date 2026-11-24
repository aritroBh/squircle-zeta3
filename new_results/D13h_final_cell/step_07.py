import math, time, json
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
# point probes at increasing N to see piece trends as y -> 0
for Nv in [400, 1e5, 1e7, 1e9]:
    yv = ys(Nv)
    t0=time.time()
    r = certified_y38(yv*(1-1e-7), yv*(1+1e-7), J=24, nsub=2)
    dt=time.time()-t0
    print(f"N={Nv:.0e}: C={float(r['C'].mid()):.2f} ({dt:.1f}s) pieces*N^1.5:",
          {k: round(float(v.mid()),3) for k,v in r.items() if k!='C'})
with open('progress.log','a') as f:
    f.write("turn2: sweep_cells.json found: 43 cells, contiguous [400,1e10], supC=617117.55. Point probes of step_38 at N=4e2..1e9 done.\n")
