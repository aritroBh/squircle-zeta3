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
dt, out, err = run('step_48.py')
print(f'--- step_48 V(3) ({dt:.1f}s) ---'); print(out); print(err[-300:] if err else '')
dt, out, err = run('step_52.py')
print(f'--- step_52 slope(4,1) ({dt:.1f}s) ---'); print(out); print(err[-300:] if err else '')
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write('step_04: rerun step_48 (V3 nonadj, rem 0) + step_52 (slope 4,1 nonadj, rem 0)\n')
