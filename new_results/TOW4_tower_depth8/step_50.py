import subprocess, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
for ds in [('7','1'), ('8','2')]:
    t0=time.time()
    r = subprocess.run([PY, work+'/prove_adj78_fast.py']+list(ds), capture_output=True, text=True, timeout=1700)
    print(ds)
    print(r.stdout.strip())
    if r.stderr: print("ERR", r.stderr[-300:])
    print('wall', round(time.time()-t0,1),'s', flush=True)
