import subprocess, time
work = '<project>/research/08_theory/TOW3_tower_depth6'
t0 = time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','prove_adj56.py','6','0'],
                   cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout.strip(), '| rc=', r.returncode, 'wall', round(time.time()-t0,1), flush=True)
