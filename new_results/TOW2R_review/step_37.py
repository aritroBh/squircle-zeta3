import os, subprocess, time
workdir = '<project>/research/08_theory/TOW2R_review'
PY = '<project>/.venv/bin/python'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
t0 = time.time()
r = subprocess.run(['nice', '-n', '10', PY, 'step_49.py'], capture_output=True, text=True,
                   timeout=1700, cwd=os.path.join(workdir,'tow2_copy'), env=env)
print(f'--- step_49 V(4) ({time.time()-t0:.1f}s) ---'); print(r.stdout); print(r.stderr[-300:] if r.stderr else '')
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_05: rerun step_49 (V4 nonadjacent) -> remainder 0\n')
