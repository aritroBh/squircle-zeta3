"""Planner check of M5b's W5'' result, independent code paths:
 (1) (3,2) Lie pairing of the symbol of gpl_w5pp.tsv (my pairing code, w4_symbol engine) and of N5'' (G-word symbols);
 (2) value: R5'' (w5pp_closed_form.json) + N5'' (cobracket5pp.json) against W5pp_1350.json, atoms by direct series."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
TH = Path(__file__).resolve().parents[2]
DPS = sys.argv[1] if len(sys.argv) > 1 else '1400'
# --- symbol part: helpers from the W6 pairing script
psrc = (TH / 'M5a_next_coefficient' / 'planner' / 'planner_w6_pairing.py').read_text()
exec(compile(psrc[:psrc.index('# weight-5 control')], 'pairing', 'exec'))
here = Path(__file__).resolve().parents[1]
S = symbol_of(load(here / 'gpl_w5pp.tsv'))
corr = json.loads((here / 'cobracket5pp.json').read_text())['correction']
def qq(pq): return Q(F(pq[0])) + Q(F(pq[1])) * R
N = {}
for t in corr:
    x, y, c = qq(t['x']), qq(t['y']), qq(t['coefficient'])
    N = add(N, scale(g_symbol(ZERO, (ZERO, ZERO, ONE / x, ZERO, ONE / (x * y)), ONE), c))
print("W5'' (3,2):", [fmt(pair(S, A3, P)) for A3 in L3])
print("N5''  (3,2):", [fmt(pair(N, A3, P)) for A3 in L3])
# --- value part: helpers from planner_eval_terms
sys.argv = ['x', str(here / 'w5pp_closed_form.json'), DPS]
esrc = (TH / 'M4e_w5_lchi' / 'planner' / 'planner_eval_terms.py').read_text()
esrc = esrc.replace('HERE = Path(__file__).resolve().parent.parent', 'HERE = TH / "M4e_w5_lchi"')
exec(compile(esrc[:esrc.index('import math\nterms =')], 'evalhelpers', 'exec'))
d = json.loads((here / 'w5pp_closed_form.json').read_text())
tot = num(F(d['aL']), F(d['bL'])) * Lchi if (F(d['aL']) or F(d['bL'])) else mp.mpf(0)
for nm, a, q in zip(d['names'], d['ca'], d['cb']):
    if F(a) or F(q):
        v = mp.mpf(1)
        for f in split(nm): v *= atom(f)
        tot += num(F(a), F(q)) * v
for t in corr:
    x = num(F(t['x'][0]), F(t['x'][1])); y = num(F(t['y'][0]), F(t['y'][1]))
    tot += num(F(t['coefficient'][0]), F(t['coefficient'][1])) * li32(x, y)
W = mp.mpf(json.loads((here / 'W5pp_1350.json').read_text())['value'])
print("aL, bL =", d['aL'], d['bL'], "; terms", sum(1 for a, q in zip(d['ca'], d['cb']) if F(a) or F(q)))
print("R5'' + N5'' - W5''(1350 file) =", mp.nstr(tot - W, 6), "at dps", DPS)
