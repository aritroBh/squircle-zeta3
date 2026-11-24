import sys, os, time, pickle
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1); os.chdir(fd1)
from fractions import Fraction as Fr
from cwfield import QSI
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain
c = QSI(0, Fr(1,2)); r = QSI(3, -2); F = QSI
red = Reducer(F, r); ci = CornerIntegrand(F, c, r)
data=pickle.load(open('kern_sq_11_20.pkl','rb'))['kernels']
for k in [21,22]:
    t0=time.time()
    ipf = build_integrand_pf(ci, red, k)
    fin = red.reduce_k(ipf)
    data[k] = {key: pf_plain(P) for key,P in fin.items()}
    print(f"squircle k={k}: {time.time()-t0:.1f}s", flush=True)
pickle.dump({'kernels':data}, open('kern_sq_11_22.pkl','wb'))
print("saved")
open(fd1+'/progress.log','a').write('step_52: squircle QSI kernels extended to k=22\n')
