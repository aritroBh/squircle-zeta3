import sys, os, time, pickle
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1); os.chdir(fd1)
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain
ANGLES = {'c23': (Fr(2,3), Fr(1,5)), 'c57': (Fr(5,7), Fr(1,6)),
          'c14': (Fr(1,4), Fr(3,5)), 'c67': (Fr(6,7), Fr(1,13))}
for name,(c,r) in ANGLES.items():
    F = make_quad(-r); red = Reducer(F, F(r)); ci = CornerIntegrand(F, F(c), F(r))
    data = pickle.load(open(f'kern20_{name}.pkl','rb'))['kernels']
    for k in [21,22]:
        t0=time.time()
        ipf = build_integrand_pf(ci, red, k)
        fin = red.reduce_k(ipf)
        data[k] = {key: pf_plain(P) for key,P in fin.items()}
        print(f"{name} k={k}: {time.time()-t0:.1f}s", flush=True)
    meta = pickle.load(open(f'kern20_{name}.pkl','rb'))
    pickle.dump({'c':meta['c'],'r':meta['r'],'kernels':data}, open(f'kern22_{name}.pkl','wb'))
    print(f"saved kern22_{name}.pkl", flush=True)
open(fd1+'/progress.log','a').write('step_61: corner kernels k=21,22 for new angles\n')
