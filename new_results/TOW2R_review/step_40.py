import os, subprocess, time
workdir = '<project>/research/08_theory/TOW2R_review'
PY = '<project>/.venv/bin/python'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
for delta in (3,4):
    t0 = time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_E3b34.py',str(delta)], capture_output=True, text=True,
                       timeout=1700, cwd=os.path.join(workdir,'tow2_copy'), env=env)
    print(f'--- prove_E3b34.py {delta} ({time.time()-t0:.1f}s) rc={r.returncode} ---')
    print(r.stdout.strip())
    if r.returncode: print(r.stderr[-400:])
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_08: rerun prove_E3b34 delta=3,4 -> True\n')
