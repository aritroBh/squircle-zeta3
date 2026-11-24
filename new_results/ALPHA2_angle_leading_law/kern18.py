"""kern18.py -- extend angles to k=17,18. args: group (old/new)"""
import sys, time, pickle, os
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain

grp = sys.argv[1]
OLD = {'c12': (Fr(1,2), Fr(1,3)), 'c35': (Fr(3,5), Fr(1,4)),
       'c13': (Fr(1,3), Fr(1,2)), 'c45': (Fr(4,5), Fr(1,9))}
NEW = {'c23': (Fr(2,3), Fr(1,5)), 'c57': (Fr(5,7), Fr(1,6)),
       'c14': (Fr(1,4), Fr(3,5)), 'c67': (Fr(6,7), Fr(1,13))}
ANGLES = OLD if grp=='old' else NEW
for name, (c, r) in ANGLES.items():
    F = make_quad(-r)
    red = Reducer(F, F(r))
    ci = CornerIntegrand(F, F(c), F(r))
    base = f'kernx_{name}.pkl' if grp=='old' else f'kern_{name}.pkl'
    data = pickle.load(open(base,'rb'))['kernels']
    for k in [17, 18]:
        t0 = time.time()
        ipf = build_integrand_pf(ci, red, k)
        fin = red.reduce_k(ipf)
        data[k] = {key: pf_plain(P) for key, P in fin.items()}
        print(f"{name} k={k}: {time.time()-t0:.1f}s", flush=True)
    pickle.dump({'c': (c.numerator, c.denominator), 'r': (r.numerator, r.denominator),
                 'kernels': data}, open(f'kern18_{name}.pkl', 'wb'))
    print(f"{name}: saved kern18_{name}.pkl", flush=True)
print("GRP DONE", grp)
