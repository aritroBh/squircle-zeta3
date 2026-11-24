import sys, os, time, pickle
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1); os.chdir(fd1)
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain
c, r = Fr(1,2), Fr(1,3)
F = make_quad(-r)
red = Reducer(F, F(r))
ci = CornerIntegrand(F, F(c), F(r))
t0=time.time()
ipf = build_integrand_pf(ci, red, 19)
fin = red.reduce_k(ipf)
print(f"k=19 build+reduce: {time.time()-t0:.1f}s, keys={len(fin)}")
