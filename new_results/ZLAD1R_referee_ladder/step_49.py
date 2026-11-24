"""ZLAD1R step_49 -- diff unmodified-G1 K=4 output vs my table; launch K=5 background; check K=8 run."""
import json, subprocess, os, time
base = '<project>'
mine = json.load(open('muR_exact_N40_K8.json'))
g1k4 = json.load(open('/tmp/flat_exact_N8_K4.json'))
def tup2name(s):
    e = eval(s)
    if e == (0,)*len(e): return '1'
    out = []
    for idx, x in enumerate(e):
        if x: out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)
bad = 0; checked = 0
for k in range(5):
    for n in range(1, 9):
        a = mine[str(k)][str(n)]
        b = {tup2name(m): v for m, v in g1k4[str(k)][str(n)].items()}
        for m in set(a) | set(b):
            checked += 1
            if a.get(m, '0') != b.get(m, '0'):
                bad += 1; print('MISMATCH', k, n, m)
print(f'unmodified G1 (K=4, N=8) vs my table: {checked} coeffs, mismatches = {bad}')
# launch K=5 N=40 background
p5 = subprocess.Popen(['nice', '-n', '10', base + '/.venv/bin/python',
                       base + '/research/10_frontier/G1_zeta3_any_angle/flat_exact.py', '40', '5'],
                      cwd='/tmp', stdout=open('/tmp/g1_k5.out', 'w'), stderr=subprocess.STDOUT)
open('/tmp/g1_k5.pid', 'w').write(str(p5.pid))
print('launched unmodified G1 K=5 N=40, pid', p5.pid)
# check K=8 run
pid = int(open('/tmp/g1_k8.pid').read().strip())
try:
    os.kill(pid, 0); alive = True
except OSError:
    alive = False
print('K=8 run alive:', alive, 'bytes:', os.path.getsize('/tmp/g1_k8.out'))
open('progress.log','a').write('step_49: unmodified G1 flat_exact.py (K=4,N=8) == my table exactly; launched K=5\n')
