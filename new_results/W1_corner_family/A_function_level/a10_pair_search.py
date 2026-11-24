"""a10_pair_search.py -- symbol-level search for a "vanishing-pair" realisation of the non-classical part.

Building block  N(x,y) = G(1/x,0,0,1/y;1),  delta22 = {x}^{y}  (a06).  Candidate terms N(X,Y) with X, Y
strictly admissible functions of c (alphabet A = {c,1-c,1+c,1+2c,1-2c}) and mixing letter X-Y in A.
For F in {U, V} (exact HPL forms, a04) look for Theta = sum_k e_k N(X_k,Y_k) with
   (i)  delta22(Theta) = delta22(F)                         (function level, exact)
   (ii) point content zero: for every ordered point pair (X(1/sqrt2), Y(1/sqrt2)) the coefficients sum to 0
        (then Theta(1/sqrt2) = 0 identically);
and test whether rho(S(F - Theta)) lies in the span of rho S(Li4(f)) over 'nice' f (point value reducible to
the conjecture basis by D4/I4/V4), using the remaining freedom ker(i) cap ker(ii).
All symbol algebra exact (Fractions / sympy rationals); point values only used as labels (30 digits).
"""
import sys, pickle, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import sympy as sp
import mpmath as mp
from sym2 import SymbolRing, delta22, add
from admissible import enumerate_admissible

mp.mp.dps = 30
c = sp.Symbol('c')
R = SymbolRing([c])
BASE = [c, 1 - c, 1 + c, 1 + 2*c, 1 - 2*c]
for f in BASE:
    R.letter_index(f)
NL = len(BASE)

# ---------- 2-variable symbol of N(x,y) = G(1/x,0,0,1/y;1) in letters x,y,1-x,1-y,x-y
x, y = sp.symbols('x y')
R2 = SymbolRing([x, y])
for f in (x, y, 1 - x, 1 - y, x - y):
    R2.letter_index(f)
SN = R2.gpl_symbol((1/x, 0, 0, 1/y))
assert len(R2.alphabet) == 5, R2.alphabet
print('N(x,y) symbol terms:', len(SN), ' 2-var alphabet:', R2.alphabet)

def vec(expr):
    """log-vector of a rational function of c over A (constants dropped); None if not in A."""
    d = R.logs(expr)
    if any(i >= NL for i in d):
        return None
    return d

def pull(S, letter_vecs):
    out = defaultdict(Fr)
    for word, v in S.items():
        for ch in itertools.product(*(letter_vecs[i].items() for i in word)):
            co = v
            for _, d in ch:
                co *= d
            out[tuple(i for i, _ in ch)] += co
    return {k: x for k, x in out.items() if x}

def rho(S):
    memo = {}
    def r(word):
        if word in memo:
            return memo[word]
        n = len(word)
        if n == 1:
            res = {word: Fr(1)}
        else:
            res = defaultdict(Fr)
            for w, v in r(word[:-1]).items():
                res[w + (word[-1],)] += Fr(n - 1, n)*v
            for w, v in r(word[1:]).items():
                res[w + (word[0],)] -= Fr(n - 1, n)*v
            res = {k: v for k, v in res.items() if v}
        memo[word] = res
        return res
    out = defaultdict(Fr)
    for w, v in S.items():
        for u, xx in r(w).items():
            out[u] += v*xx
    return {k: v for k, v in out.items() if v}

def li4_symbol(f):
    A = R.logs(1 - f); B = R.logs(f)
    out = defaultdict(Fr)
    for i, a in A.items():
        for j, b in B.items():
            for k, cc in B.items():
                for l, d in B.items():
                    out[(i, j, k, l)] += -a*b*cc*d
    return {k: v for k, v in out.items() if v}

# ---------- admissible functions
E = int(sys.argv[1]) if len(sys.argv) > 1 else 2
roots = [Fr(0), Fr(1), Fr(-1), Fr(-1, 2), Fr(1, 2)]
adm = []
for s, k, ex in enumerate_admissible(roots, E=E, K=3):
    f = sp.Integer(s)*sp.Integer(2)**k
    for r, e in zip(roots, ex):
        f *= (c - sp.Rational(r.numerator, r.denominator))**e
    adm.append(sp.cancel(f))
r2 = mp.sqrt(2)
def ptv(f):
    return sp.lambdify(c, f, 'mpmath')(1/r2)
admv = [(f, ptv(f)) for f in adm]
print('admissible functions:', len(adm), flush=True)

# ---------- candidate N-terms
terms = []
for (X, xv), (Y, yv) in itertools.product(admv, repeat=2):
    if X == Y:
        continue
    vm = vec(X - Y)
    if vm is None:
        continue
    lv = [vec(X), vec(Y), vec(1 - X), vec(1 - Y), vm]
    if any(v is None for v in lv):
        continue
    S = pull(SN, lv)
    terms.append({'X': X, 'Y': Y, 'pt': (mp.nstr(xv, 20), mp.nstr(yv, 20)), 'ob': delta22(S), 'rho': rho(S)})
print('candidate N-terms with mixing letter in A:', len(terms), flush=True)

# ---------- targets
D = pickle.load(open('a04_rfamily.pkl', 'rb'))
IDX = {0: R.letter_index(c), 1: R.letter_index(1 - c), -1: R.letter_index(1 + c)}
SGN = {0: 1, 1: -1, -1: 1}
def hpl_sym(Fd):
    out = defaultdict(Fr)
    for (w, mono), v in Fd.items():
        if mono != () or len(w) != 4:
            continue
        sg = 1
        for a in w:
            sg *= SGN[a]
        out[tuple(IDX[a] for a in reversed(w))] += v*sg
    return {k: v for k, v in out.items() if v}
TGT = {'U': hpl_sym(D['U']), 'V': hpl_sym(D['V'])}

# ---------- nice Li4 columns
NV = [1/r2, -1/r2, mp.mpf(1)/2, 1 - 1/r2, 2 - r2, r2 - 1, 1 - r2, -(r2 - 1)**2, (r2 - 1)**2, (r2 - 1)**4, mp.mpf(-1)]
def isnice(v):
    for t in (v, 1/v if v != 0 else None):
        if t is not None and any(abs(t - z) < mp.mpf(10)**-20 for z in NV):
            return True
    return False
licols = [(f, v, rho(li4_symbol(f))) for f, v in admv]
nice_li = [lc for lc in licols if isnice(lc[1]) and lc[2]]
bad_li = [lc for lc in licols if not isnice(lc[1]) and lc[2]]
print('Li4 columns: nice %d, bad %d' % (len(nice_li), len(bad_li)), flush=True)

pickle.dump({'terms': [{k: (str(v) if k in ('X', 'Y') else v) for k, v in t.items()} for t in terms],
             'nice': [(str(f), mp.nstr(v, 20), r) for f, v, r in nice_li],
             'bad': [(str(f), mp.nstr(v, 20), r) for f, v, r in bad_li],
             'TGT': TGT}, open('a10_data.pkl', 'wb'))
print('saved a10_data.pkl')
