import subprocess, time, os
PY = '<project>/.venv/bin/python'
procs = {}
for dps, N in [(520, 1750), (800, 2700)]:
    log = open(f'build_{dps}.log','w')
    p = subprocess.Popen(['nice', PY, 'step_28_build_tables.py', str(dps), str(N), '10'],
                         stdout=log, stderr=subprocess.STDOUT, cwd='.')
    procs[dps] = p
print("launched:", {k: v.pid for k, v in procs.items()})
# wait for the 520 one
t0 = time.time()
while time.time()-t0 < 1400:
    time.sleep(15)
    st = {d: p.poll() for d, p in procs.items()}
    print("%.0fs:" % (time.time()-t0), st)
    if all(s is not None for s in st.values()):
        break
for d in procs:
    print('='*40, d)
    print(open(f'build_{d}.log').read())
