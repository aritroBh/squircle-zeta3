import mpmath as mp, sys, json
from jet import *
from polys import polys
mp.mp.dps = int(sys.argv[1]); h = mp.mpf(1)/int(sys.argv[2])
out = {}
for name, U in polys().items():
    P = perimeter_jet(U, h=h)
    D = [P[0], P[1], 2*P[2], 6*P[3]]
    out[name] = [mp.nstr(x, mp.mp.dps-5) for x in D]
    print(name, *[mp.nstr(x, 50) for x in D], flush=True)
json.dump(out, open('jet_dps%d_h%s.json' % (mp.mp.dps, sys.argv[2]), 'w'), indent=1)
