"""a10c_pair_search_4letters.py -- a10b restricted to the alphabet {c,1-c,1+c,1+2c} (roots 0,1,-1,-1/2: none on the path).
exponent vectors (no sympy factoring in the inner loop).  Letters of the alphabet A are the linear factors
(c - r), r in ROOTS = (0, 1, -1, -1/2, 1/2)  <->  c, 1-c, 1+c, 1+2c, 1-2c  (constants dropped).
An admissible f is (s, k, e): f = s 2^k prod (c - r)^{e_r}.  Output: a10_data.pkl with
  terms: N(X,Y) = G(1/X,0,0,1/Y;1) with X-Y in A: symbol pulled back from the 2-variable symbol of N,
         delta22, rho, point labels (values at c = 1/sqrt2), and the exact (s,k,e) data of X, Y;
  nice / bad Li4 columns: rho S(Li4(f)).
"""
import sys, pickle, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import sympy as sp
import mpmath as mp
from sym2 import SymbolRing, delta22
from admissible import enumerate_admissible, pmul, psub, factor_over, letter_poly

mp.mp.dps = 30
ROOTS = [Fr(0), Fr(1), Fr(-1), Fr(-1, 2)]
NR = len(ROOTS)
E = int(sys.argv[1]) if len(sys.argv) > 1 else 2

def poly_of(s, k, e, part):            # numerator (part='N') or denominator (part='D') polynomial
    P = [Fr(s)*Fr(2)**k] if part == 'N' else [Fr(1)]
    for r, ex in zip(ROOTS, e):
        if (part == 'N' and ex > 0) or (part == 'D' and ex < 0):
            for _ in range(abs(ex)):
                P = pmul(P, letter_poly(r))
    return P

def dvec(e):
    return tuple(max(-x, 0) for x in e)

def mvec(mult):
    return tuple(mult.get(r, 0) for r in ROOTS)

def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))

adm = []
for s, k, e in enumerate_admissible(ROOTS, E=E, K=3):
    N, D = poly_of(s, k, e, 'N'), poly_of(s, k, e, 'D')
    res = factor_over(psub(D, N), ROOTS)
    assert res is not None
    one_minus = vsub(mvec(res[1]), dvec(e))
    adm.append({'s': s, 'k': k, 'e': e, 'N': N, 'D': D, 'v': e, 'v1': one_minus})
r2 = mp.sqrt(2); x0 = 1/r2
def val(f):
    v = mp.mpf(f['s'])*mp.mpf(2)**f['k']
    for r, ex in zip(ROOTS, f['e']):
        v *= (x0 - mp.mpf(r.numerator)/r.denominator)**ex
    return v
def val0(f):     # value at c = 0 (None if 0 or infinite)
    if f['e'][0] != 0:
        return 'zero' if f['e'][0] > 0 else 'inf'
    v = Fr(f['s'])*Fr(2)**f['k']
    for r, ex in zip(ROOTS[1:], f['e'][1:]):
        v *= (-r)**ex
    return v
for f in adm:
    f['pt'] = val(f); f['at0'] = val0(f)
print('admissible:', len(adm), flush=True)

# 2-variable symbol of N(x,y) = G(1/x,0,0,1/y;1)
x, y = sp.symbols('x y')
R2 = SymbolRing([x, y])
for g in (x, y, 1 - x, 1 - y, x - y):
    R2.letter_index(g)
SN = R2.gpl_symbol((1/x, 0, 0, 1/y))
assert len(R2.alphabet) == 5

def tovec(t):
    return {i: v for i, v in enumerate(t) if v}

def pull(S, letter_vecs):
    out = defaultdict(Fr)
    for word, v in S.items():
        for ch in itertools.product(*(letter_vecs[i].items() for i in word)):
            co = Fr(v)
            for _, d in ch:
                co *= d
            out[tuple(i for i, _ in ch)] += co
    return {k: x for k, x in out.items() if x}

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
        for u, xx in r(w).items(): out[u] += v*xx
    return {k: v for k, v in out.items() if v}

terms = []
for fX, fY in itertools.product(adm, repeat=2):
    if fX is fY:
        continue
    # X - Y = (N_X D_Y - N_Y D_X) / (D_X D_Y)
    num = psub(pmul(fX['N'], fY['D']), pmul(fY['N'], fX['D']))
    res = factor_over(num, ROOTS)
    if res is None:
        continue
    vm = vsub(vsub(mvec(res[1]), dvec(fX['e'])), dvec(fY['e']))
    lv = [tovec(fX['v']), tovec(fY['v']), tovec(fX['v1']), tovec(fY['v1']), tovec(vm)]
    S = pull(SN, lv)
    terms.append({'X': (fX['s'], fX['k'], fX['e']), 'Y': (fY['s'], fY['k'], fY['e']),
                  'pt': (mp.nstr(fX['pt'], 20), mp.nstr(fY['pt'], 20)),
                  'at0': (fX['at0'], fY['at0']),
                  'ob': delta22(S), 'rho': rho(S)})
print('N-terms with mixing in A:', len(terms), flush=True)

def li4sym(f):
    A = tovec(f['v1']); B = tovec(f['v'])
    out = defaultdict(Fr)
    for i, a in A.items():
        for j, b in B.items():
            for k, cc in B.items():
                for l, d in B.items():
                    out[(i, j, k, l)] += -a*b*cc*d
    return {k: v for k, v in out.items() if v}
NV = [1/r2, -1/r2, mp.mpf(1)/2, 1 - 1/r2, 2 - r2, r2 - 1, 1 - r2, -(r2 - 1)**2, (r2 - 1)**2, (r2 - 1)**4, mp.mpf(-1)]
def isnice(v):
    for t in (v, 1/v):
        if any(abs(t - z) < mp.mpf(10)**-20 for z in NV):
            return True
    return False
nice, bad = [], []
for f in adm:
    rr = rho(li4sym(f))
    if not rr:
        continue
    item = ((f['s'], f['k'], f['e']), mp.nstr(f['pt'], 20), rr, f['at0'])
    (nice if isnice(f['pt']) else bad).append(item)
print('Li4 columns: nice %d, bad %d' % (len(nice), len(bad)), flush=True)

# targets: U, V symbols in this letter indexing (letters: 0:c, 1:1-c, 2:1+c, 3:1+2c, 4:1-2c)
D = pickle.load(open('a04_rfamily.pkl', 'rb'))
IDX = {0: 0, 1: 1, -1: 2}
SGN = {0: 1, 1: -1, -1: 1}
def hpl_sym(Fd):
    out = defaultdict(Fr)
    for (w, mono), v in Fd.items():
        if mono != () or len(w) != 4: continue
        sg = 1
        for a in w: sg *= SGN[a]
        out[tuple(IDX[a] for a in reversed(w))] += v*sg
    return {k: v for k, v in out.items() if v}
TGT = {'U': hpl_sym(D['U']), 'V': hpl_sym(D['V'])}
pickle.dump({'terms': terms, 'nice': nice, 'bad': bad, 'TGT': TGT, 'ROOTS': ROOTS}, open('a10_data4E3.pkl', 'wb'))
print('saved a10_data.pkl', flush=True)
