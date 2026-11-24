import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import subprocess, sys
t0 = time.time()
r = subprocess.run(['nice','-n','10',sys.executable,'prove_S34.py','4'], capture_output=True, text=True, timeout=2800)
print(r.stdout.strip(), r.stderr[-150:] if r.returncode else '', f'wall {round(time.time()-t0,1)}s', flush=True)
