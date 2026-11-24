
import sys, time
sys.path.insert(0,'.')
from genintegrand import S_poly
for d in [10,12,14,16,18,20]:
    t0=time.time()
    Ss = [S_poly(d,m) for m in range(0,d+1,5)]
    print(f"d={d}: 5 S_polys in {time.time()-t0:.1f}s")
