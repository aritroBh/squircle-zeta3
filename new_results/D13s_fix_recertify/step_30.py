
import shutil, os, subprocess
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
shutil.copy(os.path.join(WORK,'D13h_final_cell__step_18.py'), os.path.join(FIX,'step_18.py'))
shutil.copy(os.path.join(WORK,'D13h_final_cell__planner_sweep_ratio.py'), os.path.join(FIX,'planner_sweep_ratio.py'))

# smoke test: two cells via the fixed certified_y38
test = '''
import math, time
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
for (n1,n2) in [(100.0,105.0),(390.0,400.0),(400.0,600.0),(1e6,1.5e6)]:
    t0=time.time()
    r = certified_y38(ys(n2), ys(n1), J=24, nsub=4)
    print(f"[{n1},{n2}] C = {float(r[chr(67)].upper()):.4f}  ({time.time()-t0:.1f}s)")
'''
open(os.path.join(FIX,'smoke.py'),'w').write(test)
r = subprocess.run([PY, 'smoke.py'], cwd=FIX, capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-1500:] if r.stderr else "")
