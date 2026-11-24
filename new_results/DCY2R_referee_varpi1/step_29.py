import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
t0=time.monotonic()
from mpmath import mp
mp.dps = 40
import dcy2_audit as A
# rerun their audit at more N, especially the worst cell region
for N in [630, 3075, 3137, 5000, 30000, 100000]:
    A.audit(N)
print("elapsed", time.monotonic()-t0)
