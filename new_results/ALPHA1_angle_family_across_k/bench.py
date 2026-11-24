import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf
c = Fr(1,2); r = Fr(1,3)
F = make_quad(-r)
red = Reducer(F, F(r))
ci = CornerIntegrand(F, F(c), F(r))
import sys as s2
kmax = int(s2.argv[1])
for k in range(9, kmax+1):
    t0 = time.time()
    ipf = build_integrand_pf(ci, red, k)
    fin = red.reduce_k(ipf)
    print(f"k={k}: {time.time()-t0:.1f}s keys={len(fin)}", flush=True)
