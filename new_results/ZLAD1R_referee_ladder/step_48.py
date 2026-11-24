import subprocess, os, time
base = '<project>'
g1 = base + '/research/10_frontier/G1_zeta3_any_angle/flat_exact.py'
t0 = time.time()
r = subprocess.run(['nice', '-n', '10', base + '/.venv/bin/python', g1, '8', '4'],
                   cwd='/tmp', capture_output=True, text=True, timeout=900)
print(f'K=4 N=8 unmodified G1 (abs path): {time.time()-t0:.0f}s, exit {r.returncode}')
print(r.stdout[-800:])
print('err:', r.stderr[-300:])
