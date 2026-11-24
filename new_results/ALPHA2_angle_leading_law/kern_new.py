"""step_08.py -- exact kernels for NEW angles c in {2/3, 5/7, 1/4, 6/7}, k = 1..16, saved to disk."""
import sys, time, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
import os
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain

ANGLES = {'c23': (Fr(2,3), Fr(1,5)), 'c57': (Fr(5,7), Fr(1,6)),
          'c14': (Fr(1,4), Fr(3,5)), 'c67': (Fr(6,7), Fr(1,13))}
KMAX = 16
for name, (c, r) in ANGLES.items():
    F = make_quad(-r)
    red = Reducer(F, F(r))
    ci = CornerIntegrand(F, F(c), F(r))
    data = {}
    t00 = time.time()
    for k in range(1, KMAX+1):
        t0 = time.time()
        ipf = build_integrand_pf(ci, red, k)
        fin = red.reduce_k(ipf)
        data[k] = {key: pf_plain(P) for key, P in fin.items()}
        print(f"{name} k={k}: {time.time()-t0:.1f}s", flush=True)
    pickle.dump({'c': (c.numerator, c.denominator), 'r': (r.numerator, r.denominator),
                 'kernels': data}, open(f'kern_{name}.pkl', 'wb'))
    print(f"{name}: total {time.time()-t00:.1f}s, saved kern_{name}.pkl", flush=True)
print("STEP08 DONE")
