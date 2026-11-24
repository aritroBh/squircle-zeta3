"""PUR2 (C),(E): for kappa in (kappa_D, kappa_max): proved tau0 and numerical sharp tau*, and where D<0."""
import mpmath as mp, sys
sys.path.insert(0, '.')
mp.mp.dps = 30
from pur2_numerics_lib import Curve, PI4
Ng = 300
grid = [PI4 * i / Ng for i in range(1, Ng + 1)]
for kap in ['1.3', '1.36', '1.4', '1.45', '1.48', '1.5', '1.507']:
    kap = mp.mpf(kap)
    C = Curve('k', kappa=kap)
    Ds = [C.D(u)[1] for u in grid]
    neg = [grid[i] for i in range(Ng) if Ds[i] < 0]
    tau0 = 1 / (1 + kap * mp.tan(kap * C.V))
    ok = lambda tau: min(C.hpp(u, tau) for u in grid) >= 0
    if ok(mp.mpf(1)):
        ts = mp.mpf(1)
    else:
        a, b = mp.mpf('1e-4'), mp.mpf(1)
        for _ in range(30):
            m = (a + b) / 2
            a, b = (m, b) if ok(m) else (a, m)
        ts = a
    print("kappa=%s kappa^2=%s V=%s phi(V)=%s minDphi/M^2=%s D<0 on u in [%s,%s] tau0(proved)=%s tau*(grid)=%s" % (
        mp.nstr(kap, 5), mp.nstr(kap**2, 5), mp.nstr(C.V, 8), mp.nstr(C.phi(C.V), 6), mp.nstr(min(Ds), 6),
        mp.nstr(neg[0], 4) if neg else '-', mp.nstr(neg[-1], 4) if neg else '-', mp.nstr(tau0, 6), mp.nstr(ts, 6)))
    sys.stdout.flush()
