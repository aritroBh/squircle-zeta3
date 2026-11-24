import subprocess, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
for ds in [('8','4'), ('7','2'), ('8','3')]:
    t0=time.time()
    r = subprocess.run([PY, work+'/prove_adj78.py']+list(ds), capture_output=True, text=True, timeout=1800)
    print(r.stdout.strip(), r.stderr[-300:] if r.stderr else '', 'wall', round(time.time()-t0,1),'s', flush=True)
