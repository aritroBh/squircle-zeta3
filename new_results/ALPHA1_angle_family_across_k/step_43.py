step08 = r'''"""step_08.py -- exact kernels for c in {1/2, 3/5, 1/3, 4/5}, k = 1..14, saved to disk."""
import sys, time, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf, pf_plain

ANGLES = {'c12': (Fr(1,2), Fr(1,3)), 'c35': (Fr(3,5), Fr(1,4)),
          'c13': (Fr(1,3), Fr(1,2)), 'c45': (Fr(4,5), Fr(1,9))}
KMAX = 14
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
'''
open('step_08.py','w').write(step08)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'step_08.py'], capture_output=True, text=True, timeout=1700)
print("elapsed", time.time()-t0)
print(r.stdout[-2500:])
print(r.stderr[-1000:] if r.returncode else "RC=0")
