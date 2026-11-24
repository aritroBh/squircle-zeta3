"""Planner: reduced representative of R5_nice (lchi_final_1200.json).
Li5(b^2) = 16 Li5(b) + 16 Li5(-b) and Li5(-1) = -(15/16) zeta(5) (exact identities); drops zero terms.
Writes planner/R5_nice_reduced.json and prints the Li5 / zeta(5) / L(chi8,5) coefficients."""
import json, math
from fractions import Fraction as F
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
d = json.loads((HERE / 'lchi_final_1200.json').read_text())
c = {nm: [F(a), F(q)] for nm, a, q in zip(d['names'], d['ca'], d['cb'])}
B2, Bp, Bm = 'Li5((3-2*sqrt2))', 'Li5((-1+1*sqrt2))', 'Li5((1-1*sqrt2))'
a, q = c.pop(B2)
for k in (Bp, Bm): c[k][0] += 16 * a; c[k][1] += 16 * q
a, q = c.pop('Li5(-1)')
c['z5'][0] += F(-15, 16) * a; c['z5'][1] += F(-15, 16) * q
c = {nm: v for nm, v in c.items() if v[0] or v[1]}
c['Lchi8_5'] = [F(d['aL']), F(d['bL'])]
for nm, (a, q) in c.items():
    if nm.startswith('Li5') or nm in ('z5', 'Lchi8_5'): print(f'{nm:24s} {a} + ({q})*sqrt2')
dens = [x.denominator for v in c.values() for x in v]; nums = [abs(x.numerator) for v in c.values() for x in v]
print('terms incl. 4 Li32:', len(c) + 4, ' Li5:', sum(nm.startswith('Li5') for nm in c),
      ' max|num|:', max(nums), ' lcm den:', math.lcm(*dens))
json.dump({'source': 'lchi_final_1200.json, reduced by planner_reduce.py',
           'N5': 'unchanged (cobracket.json)', 'terms': {nm: [str(a), str(q)] for nm, (a, q) in c.items()}},
          open(HERE / 'planner' / 'R5_nice_reduced.json', 'w'), indent=1)
