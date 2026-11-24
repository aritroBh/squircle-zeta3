
import json, time, sys
from mpmath import mp, mpf, fabs
mp.dps = 130
exec(open('av2core.py').read())
plist = eval(sys.argv[1])
try:
    tab = json.load(open('b3_hifi.json'))
except Exception:
    tab = {}
for p in plist:
    t0 = time.time()
    FPa = FPrecompC(p, '0.004', 1375)
    vc, base, Fstar = B3_prod(p, '0.008', 690, FPa)
    vf, _, _ = B3_prod(p, '0.004', 1375, FPa)
    d1 = fabs(vf-vc)
    tab[str(p)] = dict(B3=mp.nstr(vf, 125), coarse=mp.nstr(vc,125), d=mp.nstr(d1,3))
    print("p=%7s B3=%s" % (p, mp.nstr(vf, 60)))
    print("          |fine-coarse|=%s  (%.0fs)" % (mp.nstr(d1,3), time.time()-t0), flush=True)
    json.dump(tab, open('b3_hifi.json','w'), indent=1)
print("done")
