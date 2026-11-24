
import json, time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 105
exec(open('av2core.py').read())
tab = json.load(open('b3_table.json'))
print("=== validate cform/vform against existing table ===")
for p in [1, 2, 4]:
    ref = mpf(tab[str(p)]['B3'])
    t0=time.time(); vc = B3_cform(p, '0.09', 110, '0.09', 110); t1=time.time()
    vv = B3_vform(p, '0.09', 110, '0.09', 110); t2=time.time()
    print("p=%s  cform-ref=%s (%.0fs)  vform-ref=%s (%.0fs)" % (
        p, mp.nstr(fabs(vc-ref),4), t1-t0, mp.nstr(fabs(vv-ref),4), t2-t1))
