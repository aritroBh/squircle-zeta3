"""a13_solve_filtered.py -- vanishing-pair systems with the side conditions needed for a function-level proof.

Terms N(X,Y) = G(1/X,0,0,1/Y;1) kept only if, on the path c in (0, 1/sqrt2]:
  * X, Y finite (no factor (1-2c) in a denominator) and X(c) < 1, Y(c) < 1   (N real-analytic on the path);
  * X(0) = 0 or Y(0) = 0  (N -> 0 as c -> 0), and the other one tends to 0, 1/2, -1 or -infinity(for X only);
Nice Li4 columns kept only if f < 1 on the whole path (or f > 1 on the whole path), f finite on the path and
f(0) in {0, 1/2, -1, 1} (level-2 boundary constants).
Columns are ordered by simplicity so that the rref basic solution prefers simple terms.
Usage: a13_solve_filtered.py U|V [datafile]
"""
import sys, pickle
from fractions import Fraction as Fr
from collections import defaultdict
import flint, mpmath as mp
from sym2 import delta22

DATA = sys.argv[2] if len(sys.argv) > 2 else 'a10_data.pkl'
D = pickle.load(open(DATA, 'rb'))
ROOTS = D['ROOTS']
name = sys.argv[1] if len(sys.argv) > 1 else 'V'
mp.mp.dps = 30

def fval(f, x):
    s, k, e = f
    v = mp.mpf(s)*mp.mpf(2)**k
    for r, ex in zip(ROOTS, e):
        v *= (x - mp.mpf(r.numerator)/r.denominator)**ex
    return v

GRID = [mp.mpf(i)/1000 for i in range(1, 708)] + [1/mp.sqrt(2)]
def path_ok(f, allow_above=False):
    s, k, e = f
    if len(e) > 4 and e[4] < 0:       # (c - 1/2) in denominator: pole on the path
        return False
    vals = [fval(f, x) for x in GRID]
    below = all(v < 1 - mp.mpf(10)**-12 for v in vals)
    above = all(v > 1 + mp.mpf(10)**-12 for v in vals)
    return below or (allow_above and above)

def at0(f):
    s, k, e = f
    if e[0] > 0: return 0
    if e[0] < 0: return 'inf'
    v = Fr(s)*Fr(2)**k
    for r, ex in zip(ROOTS[1:], e[1:]):
        v *= (-r)**ex
    return v

OKVAL = {Fr(0), Fr(1, 2), Fr(-1)}
MZVVAL = {Fr(0), Fr(1)}
MODE = sys.argv[3] if len(sys.argv) > 3 else 'both'
def term_ok(t):
    X, Y = t['X'], t['Y']
    if MODE in ('both', 'path', 'mzv', 'lev2', 'lev2b') and not (path_ok(X) and path_ok(Y)):
        return False
    if MODE == 'path':
        return True
    if MODE == 'lev2':
        a, b = at0(X), at0(Y)
        if a == 'inf' or b == 'inf':
            return False
        S = {Fr(a), Fr(b)}
        return S <= {Fr(0), Fr(1), Fr(1, 2)} or S <= {Fr(0), Fr(1), Fr(-1)}
    if MODE == 'lev2b':
        a, b = at0(X), at0(Y)
        if a == 'inf' or b == 'inf':
            return False
        S = {Fr(a), Fr(b)}
        return S <= {Fr(0), Fr(1, 2)} or S <= {Fr(0), Fr(-1)}
    if MODE == 'mzv':
        a, b = at0(X), at0(Y)
        return a != 'inf' and b != 'inf' and Fr(a) in MZVVAL and Fr(b) in MZVVAL
    a, b = at0(X), at0(Y)
    if b == 'inf':
        return False
    if a == 'inf':
        return False                  # X -> infinity at c = 0 gives divergent sub-word boundary values
    return (a == 0 and b in OKVAL) or (b == 0 and a in OKVAL)

def cx(f):
    return sum(abs(x) for x in f[2]) + abs(f[1])

terms = [t for t in D['terms'] if term_ok(t)]
terms.sort(key=lambda t: (cx(t['X']) + cx(t['Y']), str(t['X']), str(t['Y'])))
nice = []
for f, v, rr, a0 in D['nice']:
    if MODE in ('both', 'path', 'mzv', 'lev2', 'lev2b') and not path_ok(f, allow_above=True):
        continue
    if a0 == 'inf' or (a0 != 'zero' and Fr(a0) not in (OKVAL | {Fr(1)})):
        continue
    if MODE == 'mzv' and a0 != 'zero' and Fr(a0) not in MZVVAL:
        continue
    nice.append((f, v, rr))
nice.sort(key=lambda t: cx(t[0]))
print(f'{name}: filtered N-terms {len(terms)} (of {len(D["terms"])}), filtered nice Li4 {len(nice)} (of {len(D["nice"])})', flush=True)

SF = D['TGT'][name]
def rho(S):
    memo = {}
    def r(word):
        if word in memo: return memo[word]
        n = len(word)
        if n == 1: res = {word: Fr(1)}
        else:
            res = defaultdict(Fr)
            for w, v in r(word[:-1]).items(): res[w + (word[-1],)] += Fr(n - 1, n)*v
            for w, v in r(word[1:]).items(): res[w + (word[0],)] -= Fr(n - 1, n)*v
            res = {k: v for k, v in res.items() if v}
        memo[word] = res; return res
    out = defaultdict(Fr)
    for w, v in S.items():
        for u, x in r(w).items(): out[u] += v*x
    return {k: v for k, v in out.items() if v}
obF, rhoF = delta22(SF), rho(SF)
cols = [('N', t) for t in terms] + [('L', n) for n in nice]
ob_keys = sorted(set(obF) | set().union(*[set(t['ob']) for t in terms]))
labels = sorted({t['pt'] for t in terms})
lab_idx = {q: i for i, q in enumerate(labels)}
rho_keys = sorted(set(rhoF) | set().union(*[set(t['rho']) for t in terms]) | set().union(*[set(n[2]) for n in nice]))
nrow = len(ob_keys) + len(labels) + len(rho_keys); ncol = len(cols)
M = flint.fmpq_mat(nrow, ncol + 1)
def setv(i, j, v):
    v = Fr(v); M[i, j] = flint.fmpq(v.numerator, v.denominator)
oi = {k: i for i, k in enumerate(ob_keys)}
ri = {k: len(ob_keys) + len(labels) + i for i, k in enumerate(rho_keys)}
for j, (kind, t) in enumerate(cols):
    if kind == 'N':
        for k, v in t['ob'].items(): setv(oi[k], j, v)
        setv(len(ob_keys) + lab_idx[t['pt']], j, 1)
        for k, v in t['rho'].items(): setv(ri[k], j, v)
    else:
        for k, v in t[2].items(): setv(ri[k], j, v)
for k, v in obF.items(): setv(oi[k], ncol, v)
for k, v in rhoF.items(): setv(ri[k], ncol, v)
R, rank = M.rref()
piv = []
for i in range(rank):
    for j in range(ncol + 1):
        if R[i, j] != 0:
            piv.append(j); break
ok = ncol not in piv
print('equations', nrow, 'unknowns', ncol, 'rank', rank, ' solvable:', ok, flush=True)
if ok:
    sol = [Fr(0)]*ncol
    for i, j in enumerate(piv):
        q = R[i, ncol]; sol[j] = Fr(int(q.p), int(q.q))
    used = [(j, sol[j]) for j in range(ncol) if sol[j] != 0]
    print('basic solution: %d N-terms, %d Li4 terms' % (sum(1 for j, _ in used if cols[j][0] == 'N'),
                                                        sum(1 for j, _ in used if cols[j][0] == 'L')))
    for j, v in used:
        kind, t = cols[j]
        if kind == 'N':
            print('   N %14s  X=%s Y=%s pt=%s at0=(%s,%s)' % (v, t['X'], t['Y'], t['pt'], at0(t['X']), at0(t['Y'])))
        else:
            print('   L %14s  Li4(%s) value %s' % (v, t[0], t[1]))
    pickle.dump({'cols': cols, 'sol': sol}, open('a13_sol_%s_%s_%s.pkl' % (name, DATA.replace('.pkl', ''), MODE), 'wb'))
