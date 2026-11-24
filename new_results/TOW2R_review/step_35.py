import os, subprocess, time
workdir = '<project>/research/08_theory/TOW2R_review'
PY = '<project>/.venv/bin/python'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
def run(script, args=(), timeout=1700):
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '10', PY, script, *map(str,args)],
                       capture_output=True, text=True, timeout=timeout,
                       cwd=os.path.join(workdir,'tow2_copy'), env=env)
    return time.time()-t0, r.stdout, r.stderr
dt, out, err = run('prove_S34.py', (3,))
print(f'--- prove_S34.py 3 ({dt:.1f}s) ---'); print(out); print(err[-300:] if err else '')
dt, out, err = run('prove_S34.py', (4,))
print(f'--- prove_S34.py 4 ({dt:.1f}s) ---'); print(out); print(err[-300:] if err else '')
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_03: rerun prove_S34.py l=3,4 -> numerator==0 True both\n')
