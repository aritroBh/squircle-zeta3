"""ZLAD1R step_47 -- launch unmodified G1 flat_exact.py at KMAX=5 NMAX=40 (bonus), and time a KMAX=4 NMAX=8 run."""
import subprocess, os, time
base = '<project>'
g1dir = base + '/research/10_frontier/G1_zeta3_any_angle'
# quick K=4 run in foreground to gauge timing
t0 = time.time()
r = subprocess.run(['nice', '-n', '10', base + '/.venv/bin/python', 'flat_exact.py', '8', '4'],
                   cwd='/tmp', capture_output=True, text=True, timeout=900)
print(f'K=4 N=8 unmodified G1: {time.time()-t0:.0f}s, exit {r.returncode}')
print(r.stdout[-600:])
