"""a12_solve.py -- exact linear algebra for the vanishing-pair realisation (data from a10_pair_search.py).

Unknowns: theta_k (N-terms N(X_k,Y_k)), mu_j (nice Li4 columns).  Equations (all exact over Q):
  (a) sum theta_k delta22(N_k) = delta22(F)
  (b) for every exact point pair q = (X(1/sqrt2), Y(1/sqrt2)):  sum_{k in q} theta_k = 0   [Theta(pt) = 0]
  (c) rho S(F) - sum theta_k rho S(N_k) - sum mu_j rho S(Li4(f_j)) = 0   [remainder nice mod products]
Usage: a12_solve.py U|V [allow_bad]
"""
import sys, pickle
from fractions import Fraction as Fr
from collections import defaultdict
import flint
from sym2 import delta22

D = pickle.load(open('a10_data.pkl', 'rb'))
name = sys.argv[1] if len(sys.argv) > 1 else 'V'
allow_bad = len(sys.argv) > 2 and sys.argv[2] == 'bad'
terms, nice, bad, TGT = D['terms'], D['nice'], D['bad'], D['TGT']
SF = TGT[name]

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
cols = [('N', t) for t in terms] + [('L', n) for n in nice] + ([('B', b) for b in bad] if allow_bad else [])
ob_keys = sorted(set(obF) | set().union(*[set(t['ob']) for t in terms]))
labels = sorted({t['pt'] for t in terms})
lab_idx = {q: i for i, q in enumerate(labels)}
rho_keys = sorted(set(rhoF) | set().union(*[set(t['rho']) for t in terms]) | set().union(*[set(n[2]) for n in nice])
                  | (set().union(*[set(b[2]) for b in bad]) if allow_bad else set()))
nrow = len(ob_keys) + len(labels) + len(rho_keys)
ncol = len(cols)
print(f'{name}: unknowns {ncol} (N-terms {len(terms)}, nice Li4 {len(nice)}{", bad Li4 %d" % len(bad) if allow_bad else ""}),'
      f' equations {nrow} (ob {len(ob_keys)}, point labels {len(labels)}, rho {len(rho_keys)})', flush=True)
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
# consistency: a pivot in the last column means no solution
piv = []
for i in range(rank):
    for j in range(ncol + 1):
        if R[i, j] != 0:
            piv.append(j); break
ok = ncol not in piv
print('rank', rank, ' solvable:', ok, flush=True)
if ok:
    sol = [Fr(0)]*ncol
    for i, j in enumerate(piv):
        q = R[i, ncol]; sol[j] = Fr(int(q.p), int(q.q))
    used = [(j, sol[j]) for j in range(ncol) if sol[j] != 0]
    print('basic solution with %d nonzero unknowns:' % len(used))
    for j, v in used:
        kind, t = cols[j]
        if kind == 'N':
            print('   N  %10s  X=%s  Y=%s   pt=%s' % (v, t['X'], t['Y'], t['pt']))
        else:
            print('   %s  %10s  Li4(%s)  [value %s]' % (kind, v, t[0], t[1]))
    pickle.dump({'cols': cols, 'sol': sol}, open('a12_sol_%s%s.pkl' % (name, '_bad' if allow_bad else ''), 'wb'))
