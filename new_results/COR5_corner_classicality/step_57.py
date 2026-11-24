import sys, os
g1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
sys.path.insert(0, g1)
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
os.chdir(g1)
import mpmath as mp
from hpl_series import H_series, Li31, Li22
mp.mp.dps = 60
x = mp.mpf(1)/2
checks = [
    ((0,1), mp.polylog(2, x)),
    ((0,-1), -mp.polylog(2, -x)),
    ((0,0,1), mp.polylog(3, x)),
    ((0,0,-1), -mp.polylog(3, -x)),
    ((0,0,0,1), mp.polylog(4, x)),
    ((0,0,0,-1), -mp.polylog(4, -x)),
    ((0,0,1,1), Li31(x, mp.mpf(1))),
    ((0,1,0,1), Li22(x, mp.mpf(1))),
    ((0,0,1,-1), -Li31(x, -mp.mpf(1))),
    ((0,1,0,-1), -Li22(x, -mp.mpf(1))),
    ((0,0,-1,1), -Li31(-x, -mp.mpf(1))),
    ((0,-1,0,1), -Li22(-x, -mp.mpf(1))),
]
for w, v in checks:
    vs = H_series(w, x)
    print(w, "diff:", mp.nstr(abs(vs - v), 3))
