"""a14_sparsify.py -- greedy minimal-support vanishing-pair solution (exact), U and V jointly, same candidate
set (filter 'lev2': N-terms real-analytic on c in (0,1/sqrt2], boundary values at c=0 in {0,1,1/2} or {0,1,-1};
nice Li4 columns real-analytic on the path with f(0) in {0,1/2,-1,1}).
Greedy: try to drop each column of the current joint support (most complicated first); keep it dropped if
BOTH systems stay solvable.  Output a14_sparse.pkl."""
import sys, pickle
from fractions import Fraction as Fr
from vp_common import term_ok, nice_ok, at0, cx, solve

DATA = sys.argv[1] if len(sys.argv) > 1 else 'a10_data4.pkl'
D = pickle.load(open(DATA, 'rb'))
ROOTS = D['ROOTS']
terms = [t for t in D['terms'] if term_ok(t, ROOTS, 'lev2')]
nice = [(f, v, rr) for (f, v, rr, a0) in D['nice'] if nice_ok((f, v, rr, a0), ROOTS)]
print('candidates: N-terms %d, nice Li4 %d' % (len(terms), len(nice)), flush=True)

def bscore(t):
    s = 0
    for v in (at0(t['X'], ROOTS), at0(t['Y'], ROOTS)):
        s += {Fr(0): 0, Fr(1): 3, Fr(1, 2): 2, Fr(-1): 2}.get(v, 5)
    return s
def score(c):
    return 0 if c[0] == 'L' else 1 + bscore(c[1]) + cx(c[1]['X']) + cx(c[1]['Y'])
cols = [('N', t) for t in terms] + [('L', n) for n in nice]
cols.sort(key=score)
SU, SV = D['TGT']['U'], D['TGT']['V']
sU, sV = solve(cols, SU), solve(cols, SV)
assert sU is not None and sV is not None
keep = [c for c, a, b in zip(cols, sU, sV) if a != 0 or b != 0]
print('initial joint support:', len(keep), flush=True)
changed = True
while changed:
    changed = False
    for c in sorted(keep, key=score, reverse=True):
        if c[0] != 'N':
            continue
        trial = [d for d in keep if d is not c]
        if solve(trial, SU) is not None and solve(trial, SV) is not None:
            keep = trial
            changed = True
    print('  support now', len(keep), flush=True)
sU, sV = solve(keep, SU), solve(keep, SV)
print('final joint support: %d N-terms, %d Li4' % (sum(1 for c in keep if c[0] == 'N'), sum(1 for c in keep if c[0] == 'L')))
for c, a, b in zip(keep, sU, sV):
    if c[0] == 'N':
        t = c[1]
        print('  N  U:%10s  V:%10s   X=%s Y=%s  pt=%s  at0=(%s,%s)' % (a, b, t['X'], t['Y'], t['pt'],
              at0(t['X'], ROOTS), at0(t['Y'], ROOTS)))
    else:
        print('  L  U:%10s  V:%10s   Li4(%s) value %s' % (a, b, c[1][0], c[1][1]))
pickle.dump({'final': keep, 'sU': sU, 'sV': sV, 'ROOTS': ROOTS}, open('a14_sparse.pkl', 'wb'))
