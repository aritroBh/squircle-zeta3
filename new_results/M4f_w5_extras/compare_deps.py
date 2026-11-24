"""M4f: compare the 14 dependency vectors found at 1200 vs 1500 digits (same Q-span? identical?).
Also compares the shortest target relations."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
import flint

def deps(path):
    d = json.loads((HERE/path).read_text())
    out = []
    for dep in d['deps']:
        assert dep['vec'][0] == 0 and len(dep['kept']) == 86
        out.append([F(c) for c in dep['vec'][1:]])
    return out, d

d12, j12 = deps('overcomplete_1200.json')
d15, j15 = deps('overcomplete_1500.json')
print(f'deps: {len(d12)} @1200, {len(d15)} @1500')
# identical as integer vectors?
same = d12 == d15
print('dep vectors IDENTICAL:', same)
if not same:
    # compare Q-spans via rref rank
    M = [[str(c) for c in v] for v in d12+d15]
    _, rk = flint.fmpq_mat(M).rref()
    print('combined rank:', rk, '(== 14 iff same span)')
# target relations
t12 = j12['target_rels'][0]; t15 = j15['target_rels'][0]
print('target rel identical:', t12['vec'] == t15['vec'], 'c0:', t12['vec'][0], t15['vec'][0])
print('resids:', t12['resid'], t15['resid'])
