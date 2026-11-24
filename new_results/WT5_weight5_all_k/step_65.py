base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step J: E3 direct Gamma-world LLL test at two precisions.
Basis: {lam^5, pi2 lam^3, pi4 lam, z3 lam^2, z3 pi2, z5}. Expect to recover the rel0 form."""
import json, math, sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
from m5b_common import lll_relations, log

def run(digits):
    mp.mp.dps = digits + 120
    lam = mp.log(2); pi2 = mp.pi**2; z3 = mp.zeta(3); z5 = mp.zeta(5)
    vals = [mp.mpf(json.loads((HERE/'e3_5800.json').read_text())['value_full']),
            lam**5, pi2*lam**3, (pi2**2)*lam, z3*lam**2, z3*pi2, z5]
    rels = lll_relations(vals, digits, f'E3 Gamma-world basis @{digits}')
    rels = [r for r in rels if r['vec'][0] != 0]
    print(f'@{digits}: {len(rels)} relations with target')
    for r in rels:
        print('  vec:', r['vec'], 'resid:', r['residual'], 'maxcoef:', r['maxcoef'])
    return rels

r1 = run(1100)
r2 = run(1300)
def canon(v): return tuple(v if v[0] > 0 else [-x for x in v])
c1 = {canon(r['vec']): r['residual'] for r in r1}
c2 = {canon(r['vec']): r['residual'] for r in r2}
common = sorted(set(c1) & set(c2))
print('identical at both precisions:', len(common))
for v in common:
    print('  IDENTICAL:', list(v), 'resids', c1[v], c2[v])
    # express: v[0]*E3 + sum v[i]*basis[i] = 0
    names = ['E3', 'lam^5', 'pi2*lam^3', 'pi4*lam', 'z3*lam^2', 'z3*pi2', 'z5']
    from fractions import Fraction as F
    c0 = v[0]
    expr = ' = '.join([f'E3'] + [f'({F(-c)/F(c0)})*{nm}' for c, nm in zip(v[1:], names[1:]) if c])
    print('  closed form:', expr)
out = dict(common=[list(v) for v in common], resids={str(list(v)): [c1[v], c2[v]] for v in common})
(HERE/'e3_gamma_lll.json').write_text(json.dumps(out, indent=1)+'\\n')
with open(HERE/'progress.log','a') as f:
    f.write(f'step_J gamma LLL: identical relations at 1100/1300: {len(common)}\\n')
'''
open(W+'wt5_J_gamma.py','w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_J_gamma.py'], cwd=W, capture_output=True, text=True, timeout=1200)
print(r.stdout); print('STDERR:', r.stderr[-400:])
