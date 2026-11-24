import os, subprocess, time
dst = '<project>/research/08_theory/M4i_w5_noL'
PY='<project>/.venv/bin/python'
with open(dst+'/progress.log','a') as fh:
    fh.write('step_22(nol_optimize2 750s): 295 restarts, no improvement; final L-free score (2487169, 301916160, 69). Accept.\n')

script = r'''"""M4i step D: export R5_noL.json (R5_best.json format) from state_noL.pkl.
Reductions: Li5(b^2)->16[Li5(b)+Li5(-b)]; Li5(-1)->-(15/16)z5;
            Li5(-1/sqrt2)-> Li5(1/2)/16 - Li5(1/sqrt2). Asserts no Lchi8_5 term.
"""
import json, pickle, sys
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nice_search as NS

sol = NS.solve(list(NS.S0), 'expnol')
names = sol['names']
NATOM = 146
i_b2 = names.index('Li5((3-2*sqrt2))'); i_bp = names.index('Li5((-1+1*sqrt2))')
i_bm = names.index('Li5((1-1*sqrt2))'); i_m1 = names.index('Li5(-1)')
i_z5 = names.index('z5'); i_r2m = names.index('Li5(-1/2*sqrt2)')
i_r2p = names.index('Li5(1/2*sqrt2)'); i_half = names.index('Li5(1/2)')
pool = json.loads((HERE/'pool.json').read_text())['pool']
s0 = set(NS.S0)
extras = [NS.W.Q(F(e['p']), F(e['q'])) for e in pool if NS.W.Q(F(e['p']), F(e['q'])) not in s0]
assert len(extras) == 7
def reduce3(v):
    v = list(v)
    for off in (0, 65):
        a = v[i_r2m+off]
        if a: v[i_half+off] += F(a, 16); v[i_r2p+off] -= a; v[i_r2m+off] = F(0)
        a = v[i_b2+off]
        if a: v[i_bp+off] += 16*a; v[i_bm+off] += 16*a; v[i_b2+off] = F(0)
        a = v[i_m1+off]
        if a: v[i_z5+off] += F(-15, 16)*a; v[i_m1+off] = F(0)
    return v
v = pickle.loads((HERE/'state_noL.pkl').read_bytes())['v']
assert v[130] == 0 and v[131] == 0
v = reduce3(v)
terms = {}
for i, nm in enumerate(names):
    a, b = v[i], v[65+i]
    if a or b: terms[nm] = [str(a), str(b)]
for k in range(7):
    a, b = v[132+k], v[139+k]
    if a or b: terms[f'Li5({extras[k]})'] = [str(a), str(b)]
assert not (v[130] or v[131]), 'L present!'
dens = [F(x).denominator for pq in terms.values() for x in pq]
nums = [abs(F(x).numerator) for pq in terms.values() for x in pq]
li5 = [nm for nm in terms if nm.startswith('Li5')]
print(f'R5_noL.json: terms={len(terms)} max|num|={max(nums)} lcm={lcm(*dens)}')
print(f'   Li5 args ({len(li5)}): {li5}')
json.dump({'source': 'M4i state_noL (L-free 12-dim affine optimum, score A, 3 reductions)',
           'N5': 'unchanged (cobracket.json)', 'terms': terms},
          open(HERE/'R5_noL.json', 'w'), indent=1)
print('saved R5_noL.json')
'''
open(dst+'/nol_export.py','w').write(script)
r = subprocess.run(['nice','-n','10',PY,'nol_export.py'], capture_output=True, text=True, cwd=dst, timeout=600)
print('rc=%d' % r.returncode); print(r.stdout); print('STDERR:', r.stderr[-1500:])
