import subprocess, time
PY = '<project>/.venv/bin/python'
for dps, N in [(1150, 3900), (1550, 5250)]:
    t0=time.time()
    r = subprocess.run(['nice', PY, 'step_30_wparts.py', str(dps), str(N)], capture_output=True, text=True, cwd='.')
    print("="*40, dps, "(%.0fs)" % (time.time()-t0))
    print('\n'.join(r.stdout.splitlines()[:9]))
    if r.stderr: print("STDERR:", r.stderr[-400:])
