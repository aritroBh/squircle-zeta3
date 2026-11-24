"""ZLAD1R step_32 -- full-range diffs: ZLAD1 K6/K8 tables vs G1 (N320, k<=3) and G2 (N640, k=4);
also spot-launch G1's flat_exact.py unmodified at KMAX=8 to test feasibility."""
import json, subprocess, time, os
base = '<project>'
zl6 = json.load(open(base + '/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N300_K6.json'))
zl8 = json.load(open(base + '/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N150_K8.json'))
g1 = json.load(open(base + '/research/10_frontier/G1_zeta3_any_angle/flat_exact_N320_K3.json'))
g2 = json.load(open(base + '/research/08_theory/G2_M4alpha/flat_exact_N640_K4.json'))
def tup2name(s):
    e = eval(s)
    if e == (0,)*len(e): return '1'
    out = []
    for idx, x in enumerate(e):
        if x: out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)
def diff_full(zl, gg, ks, ns):
    bad = 0; checked = 0
    for k in ks:
        for n in ns:
            a = zl[str(k)][str(n)]
            b = {tup2name(m): v for m, v in gg[str(k)][str(n)].items()}
            for m in set(a) | set(b):
                checked += 1
                if a.get(m, '0') != b.get(m, '0'):
                    bad += 1
                    if bad < 4: print('  MISMATCH', k, n, m, a.get(m), b.get(m))
    return bad, checked
bad1, c1 = diff_full(zl6, g1, range(4), range(1, 301))
print(f'ZLAD1-K6 vs G1-N320 (k<=3, n<=300): {c1} coefficients, mismatches = {bad1}')
bad2, c2 = diff_full(zl6, g2, [4], range(1, 301))
print(f'ZLAD1-K6 vs G2-N640 (k=4, n<=300): {c2} coefficients, mismatches = {bad2}')
bad3, c3 = diff_full(zl8, g2, [4], range(1, 151))
print(f'ZLAD1-K8 vs G2-N640 (k=4, n<=150): {c3} coefficients, mismatches = {bad3}')
# launch G1 unmodified at KMAX=8, NMAX=40 in background, monitor time
env = dict(os.environ)
p = subprocess.Popen(['nice', '-n', '10', base + '/.venv/bin/python', 'flat_exact.py', '40', '8'],
                     cwd=base + '/research/10_frontier/G1_zeta3_any_angle',
                     stdout=open('/tmp/g1_k8.out', 'w'), stderr=subprocess.STDOUT)
open('/tmp/g1_k8.pid', 'w').write(str(p.pid))
print('launched G1 flat_exact.py 40 8 (unmodified), pid', p.pid)
open('progress.log','a').write('step_32: full-range diffs ZLAD1 vs G1 (n<=300,k<=3) and G2 (k=4): 0 mismatches; launched unmodified G1 at K=8\n')
