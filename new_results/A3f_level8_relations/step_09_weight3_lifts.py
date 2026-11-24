"""step_09_weight3_lifts.py: weight-3 quotient with the new weight-2 relation lifted in."""
import sys, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from relations import words, make_relations
from flint import fmpq_mat, nmod_mat

a1 = pickle.load(open('audit_weight1.pickle', 'rb'))
a2 = pickle.load(open('audit_weight2.pickle', 'rb'))
new2 = pickle.load(open('weight2_new_relation.pickle', 'rb'))['row']

lower = {1: a1['selected'], 2: a2['selected'] + [new2]}
cols = words(3)
rows, counts = make_relations(3, lower)
print('rows:', len(rows), 'families:', counts, flush=True)

ids = {w: i for i, w in enumerate(cols)}
# mod two primes first
for p in (1000000009, 1000000021):
    M = nmod_mat(len(rows), len(cols), p)
    for i, row in enumerate(rows):
        for w, c in row:
            M[i, ids[w]] = c % p
    _, r = M.rref(inplace=True)
    print(f'mod {p}: rank={r} quotient={len(cols)-r}', flush=True)
# exact
t0 = time.monotonic()
M = fmpq_mat(len(rows), len(cols))
for i, row in enumerate(rows):
    for w, c in row:
        M[i, ids[w]] = int(c)
_, rq = M.rref(inplace=True)
print(f'exact Q: rank={rq} quotient={len(cols)-rq} (DG 27)  ({time.monotonic()-t0:.1f}s)', flush=True)
with open('weight3_with_lifts.pickle', 'wb') as f:
    pickle.dump({'cols': cols, 'rows': rows, 'rank_Q': rq, 'quotient': len(cols) - rq}, f)
