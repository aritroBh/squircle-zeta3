import sys, subprocess, os
print("tool python:", sys.executable)
PY = '<project>/.venv/bin/python'
print(subprocess.run([PY,'-c','import sys, sympy; print(sys.version.split()[0], sympy.__version__)'], capture_output=True, text=True).stdout)
# quick test: tow2lib L(7) construction time (foreground, cancel at end)
import time
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import sympy as sp
import tow2lib as T
t0=time.time()
L7 = T.L(7)
print("L(7) constructed", round(time.time()-t0,1), "s; now cancelling...")
t0=time.time()
L7c = sp.cancel(L7)
print("cancelled", round(time.time()-t0,1), "s; len", len(str(L7c)))
import pickle
pickle.dump(str(L7c), open('<project>/research/08_theory/TOW4_tower_depth8/L7_tow2lib.pkl','wb'))
