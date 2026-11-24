
import json, time, sys
from mpmath import mp, mpf, fabs
mp.dps = 125
exec(open('av2core.py').read())
plist = eval(sys.argv[1])
try:
    tab = json.load(open('b3_final.json'))
except Exception:
    tab = {}
for p in plist:
    t0 = time.time()
    big = mpf(p) >= 6
    gA = ('0.0055', 950) if big else ('0.008', 700)
    gB = ('0.008', 650) if big else ('0.011', 500)
    FPa = FPrecompC(p, *gA)
    vc, base, Fstar = B3_prod(p, '0.011', 500, FPa)
    vf, _, _ = B3_prod(p, '0.0055', 1000, FPa)
    FPb = FPrecompC(p, *gB)
    vb, _, _ = B3_prod(p, '0.011', 500, FPb)
    d1 = fabs(vf-vc); d2 = fabs(vf-vb)
    tab[str(p)] = dict(B3=mp.nstr(vf, 120), coarse=mp.nstr(vc,120), alt=mp.nstr(vb,120),
                       d_fine_coarse=mp.nstr(d1,3), d_grid=mp.nstr(d2,3))
    print("p=%7s B3=%s" % (p, mp.nstr(vf, 50)))
    print("          |fine-coarse|=%s |gridA-gridB|=%s  (%.0fs)" % (mp.nstr(d1,3), mp.nstr(d2,3), time.time()-t0), flush=True)
    json.dump(tab, open('b3_final.json','w'), indent=1)
print("done")
