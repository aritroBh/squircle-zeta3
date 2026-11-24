
import json, time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 115
exec(open('av2core.py').read())
plist = [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 512, 1024, 4096]
tab = {}
for p in plist:
    t0 = time.time()
    q2 = B2_quad2(p); q3 = B2_quad3(p)
    if p <= 96:
        s, nn = B2_series(p)
        err = max(fabs(s-q2), fabs(s-q3))
        tab[str(p)] = dict(B2=mp.nstr(s, 110), err=mp.nstr(err, 3), methods="series+quad2+quad3", terms=nn)
    else:
        err = fabs(q2-q3)
        tab[str(p)] = dict(B2=mp.nstr(q2, 110), err=mp.nstr(err, 3), methods="quad2+quad3")
    print("p=%6s  B2=%s" % (p, mp.nstr(mpf(tab[str(p)]['B2']), 45)))
    print("          err=%s  (%.1fs)" % (tab[str(p)]['err'], time.time()-t0))
json.dump(tab, open('b2_table.json','w'), indent=1)
print("saved b2_table.json")
