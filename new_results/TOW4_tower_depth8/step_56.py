import subprocess, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
for dl in ['7','8']:
    t0=time.time()
    r = subprocess.run([PY, work+'/prove_E3b78.py',dl], capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip()); print(r.stderr[-200:] if r.stderr else '', 'wall', round(time.time()-t0,1), flush=True)
