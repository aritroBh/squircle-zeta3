script = '''"""kern_old_ext.py -- extend old angles c12,c35,c13,c45 to k=15,16 and merge into kern_*x.pkl"""
import sys, time, pickle, os
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain

ANGLES = {'c12': (Fr(1,2), Fr(1,3)), 'c35': (Fr(3,5), Fr(1,4)),
          'c13': (Fr(1,3), Fr(1,2)), 'c45': (Fr(4,5), Fr(1,9))}
for name, (c, r) in ANGLES.items():
    F = make_quad(-r)
    red = Reducer(F, F(r))
    ci = CornerIntegrand(F, F(c), F(r))
    data = pickle.load(open(f'kern_{name}.pkl','rb'))['kernels']
    for k in [15, 16]:
        t0 = time.time()
        ipf = build_integrand_pf(ci, red, k)
        fin = red.reduce_k(ipf)
        data[k] = {key: pf_plain(P) for key, P in fin.items()}
        print(f"{name} k={k}: {time.time()-t0:.1f}s", flush=True)
    pickle.dump({'c': (c.numerator, c.denominator), 'r': (r.numerator, r.denominator),
                 'kernels': data}, open(f'kernx_{name}.pkl', 'wb'))
    print(f"{name}: saved kernx_{name}.pkl", flush=True)
print("EXT DONE")
'''
A2='<project>/research/08_theory/ALPHA2_angle_leading_law'
open(f'{A2}/kern_old_ext.py','w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
out = open(f'{A2}/kern_old_ext.out','w')
proc = subprocess.Popen(['nice','-n','10',PY,f'{A2}/kern_old_ext.py'], stdout=out, stderr=subprocess.STDOUT)
print("launched", proc.pid)
with open(f'{A2}/progress.log','a') as fh: fh.write('step: kernels 4 new angles k<=16 done (~250s each); old-angle ext launched\n')
