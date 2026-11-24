
import json, time, sys
from mpmath import mp, mpf, zeta, fabs
mp.dps = 105
exec(open('av2core.py').read())
which = sys.argv[1] if len(sys.argv) > 1 else '1'
if which == '1':
    plist = [1, 1.25, 1.5, 2, 2.5, 3, 4, 6]
else:
    plist = [8, 12, 16, 24, 32, 64, 128, 256, 1024, 4096]
try:
    tab = json.load(open('b3_table.json'))
except Exception:
    tab = {}
for p in plist:
    t0 = time.time()
    vB1 = B3_tsB(p, '0.025', 190)
    vB2 = B3_tsB(p, '0.035', 130)
    vA  = B3_tsA(p, '0.07', 130)
    errB = fabs(vB1 - vB2); errA = fabs(vB1 - vA)
    tab[str(p)] = dict(B3=mp.nstr(vB1, 100), errB=mp.nstr(errB, 3), errA=mp.nstr(errA, 3))
    print("p=%6s B3=%s" % (p, mp.nstr(vB1, 42)))
    print("         |B(h1)-B(h2)|=%s  |B-A|=%s  (%.1fs)" % (mp.nstr(errB,3), mp.nstr(errA,3), time.time()-t0), flush=True)
    json.dump(tab, open('b3_table.json','w'), indent=1)
print("saved")
