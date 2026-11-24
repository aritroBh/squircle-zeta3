import os, subprocess, time
workdir = '<project>/research/08_theory/TOW2R_review'
PY = '<project>/.venv/bin/python'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
for (delta, s) in [(4,1),(4,0)]:
    t0 = time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_adj34.py',str(delta),str(s)], capture_output=True, text=True,
                       timeout=1700, cwd=os.path.join(workdir,'tow2_copy'), env=env)
    print(f'--- prove_adj34.py {delta} {s} ({time.time()-t0:.1f}s) rc={r.returncode} ---')
    print(r.stdout.strip())
    if r.returncode: print(r.stderr[-400:])
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_07: rerun prove_adj34 (4,1),(4,0) -> True\n')
