"""Planner check of M4i's L_relation.json: L(chi8,5) = sum (p + q sqrt2) atom, with no L on the right-hand side.
Independent of the M4c/M4e/M4f/M4i code: atoms by direct series (planner_eval_terms helpers), L(chi8,5) by Hurwitz zeta.
Usage: python planner_check_Lrel.py <L_relation.json> [DPS]"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '1700']
src = Path(__file__).with_name('planner_eval_terms.py').read_text()
src = src[:src.index('import math\nterms =')]          # reuse helpers only (atoms, Lchi, num, split)
exec(compile(src.replace("TERMS = Path(sys.argv[1]).resolve()", "TERMS = Path(sys.argv[1]).resolve()"), 'helpers', 'exec'))
terms = json.loads(TERMS.read_text())['terms']
assert 'Lchi8_5' not in terms
rhs = mp.mpf(0)
for nm, (p, q) in terms.items():
    v = mp.mpf(1)
    for f in split(nm): v *= atom(f)
    rhs += num(F(p), F(q)) * v
pure_sqrt2 = all(F(p) == 0 for p, q in terms.values())
dens = [F(x).denominator for v in terms.values() for x in v]
print(f'{TERMS.name}: {len(terms)} terms, all coefficients pure sqrt2: {pure_sqrt2}, '
      f'max |num| {max(abs(F(x).numerator) for v in terms.values() for x in v)}, lcm den {math.lcm(*dens)}')
print('RHS - L(chi8,5)[Hurwitz] =', mp.nstr(rhs - Lchi, 6))
