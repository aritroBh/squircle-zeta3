code = r'''"""step_08_weight2_certificate.py: find minimal five-term witnesses; verify LLL hint is spanned;
save exact weight-2 certificate."""
import sys, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from fractions import Fraction as F
from relations import words
from flint import fmpq_mat

cols = words(2)
ids = {w: i for i, w in enumerate(cols)}
known = pickle.load(open('audit_weight2.pickle', 'rb'))['rows']
extra = pickle.load(open('weight2_fiveterm.pickle', 'rb'))
rows, names = extra['rows'], extra['names']

def mat(rowlist):
    M = fmpq_mat(len(rowlist), len(cols))
    for i, row in enumerate(rowlist):
        for w, c in row:
            M[i, ids[w]] = int(c)
    return M

M0 = mat(known); _, r0 = M0.rref(inplace=True)
singles = []
for row, nm in zip(rows, names):
    M = mat(known + [row]); _, r = M.rref(inplace=True)
    if r > r0:
        singles.append((nm, row))
        print(f'NEW RANK from: {nm}')
print(f'{len(singles)} single rows add rank (rank {r0}->{r0+1})', flush=True)

# choose the first five-term single as the certified new relation
new_row, new_name = None, None
for nm, row in singles:
    if nm.startswith('Rogers5T'):
        new_row, new_name = row, nm
        break
if new_row is None:
    new_row, new_name = singles[0][1], singles[0][0]
print('chosen new relation:', new_name)
print('terms:')
for w, c in new_row:
    print(f'   {c:+d}  H{w}')

# verify LLL hint vector is in the span of known + new_row
lll = {(6, 7): -6, (7, 0): -1, (7, 1): 7, (7, 2): 3, (7, 3): -9,
       (7, 4): -3, (7, 5): 21, (7, 6): -5, (7, 7): 5}
M = mat(known + [new_row]); rr, r = M.rref(inplace=True)
# membership test: add lll row, rank must stay r
g = F(0)
from relations import primitive
lll_row = primitive({w: F(c) for w, c in lll.items()})
M2 = mat(known + [new_row, lll_row]); _, r2 = M2.rref(inplace=True)
print(f'rank known+new={r}; rank known+new+LLLhint={r2} (equal => LLL hint spanned)', flush=True)
assert r2 == r

# quotient now 9?
print(f'quotient with new relation: {len(cols)-r} (DG benchmark 9)')
assert len(cols) - r == 9

# Solve: express LLL hint as Q-combination of (known + new_row)? Instead, relate new_row to
# LLL hint: both span the 1-dim kernel in the quotient; find scalar: reduce both to normal forms.
pp = []; start = 0
for i in range(r):
    while not rr[i, start]: start += 1
    pp.append(start); start += 1
free = [j for j in range(len(cols)) if j not in set(pp)]
def normal_form(row):
    v = {w: F(c) for w, c in row}
    nf = []
    for j in free:
        nf.append(v.get(cols[j], F(0)) - sum((v.get(cols[q], F(0)) * F(rr[i, j]) for i, q in enumerate(pp)), F(0)))
    return nf
nf_new = normal_form(new_row)
nf_lll = normal_form(lll_row)
k = next(i for i, x in enumerate(nf_new) if x != 0)
lam = nf_lll[k] / nf_new[k]
print('LLL hint =', lam, 'x new_row (in the quotient); consistent:',
      all(nf_lll[i] == lam * nf_new[i] for i in range(len(free))))

with open('weight2_new_relation.pickle', 'wb') as f:
    pickle.dump({'name': new_name, 'row': new_row, 'all_new_rows': [r for r, n in singles],
                 'all_names': [n for n, r in singles]}, f)
json.dump({'name': new_name, 'row': [[list(w), int(c)] for w, c in new_row]},
          open('weight2_new_relation.json', 'w'), indent=1)
print('saved weight2_new_relation.{pickle,json}')
'''
open('step_08_weight2_certificate.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_08_weight2_certificate.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout); print('STDERR:', r.stderr[-3000:])
