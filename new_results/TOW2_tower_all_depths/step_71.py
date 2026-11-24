import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import subprocess, sys, time
for (delta, s) in [(3,0),(4,1)]:
    t0 = time.time()
    r = subprocess.run(['nice','-n','10',sys.executable,'prove_adj34.py',str(delta),str(s)], capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip(), r.stderr[-150:] if r.returncode else '', f'wall {round(time.time()-t0,1)}s', flush=True)
