
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 100
exec(open('av2core.py').read())
ref1 = mpf('0.246617302384923740565168068839592055997186102')  # mpmath quad, dps=70, A&B agreed
print("=== p=1, fine h scan, compare vs ref ===")
for (h, n) in [(0.08, 130), (0.06, 170), (0.05, 200), (0.04, 250)]:
    t0=time.time(); v = B3_tsB(1, h, n); t1=time.time()
    print("h=%.3f n=%3d  err=%s (%.1fs)" % (h, n, mp.nstr(fabs(v-ref1),4), t1-t0))
    print("   ", mp.nstr(v, 60))
