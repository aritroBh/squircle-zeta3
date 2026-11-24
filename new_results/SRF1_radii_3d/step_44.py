
import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/claude')
import mpmath as mp
from superlib import S_area, Qfun, F2, H, h

s3 = '<project>/research/07_push/S3_superellipsoid/claude'
sv = json.load(open(s3+'/S_values.json'))

mp.mp.dps = 80
t0 = time.time()
v10 = S_area(10)
t1 = time.time()
print(f"S(10) computed in {t1-t0:.1f}s:", mp.nstr(v10, 75))
print("precomputed:          ", sv['10'])
print("match:", mp.nstr(abs(v10 - mp.mpf(sv['10'])), 3))
