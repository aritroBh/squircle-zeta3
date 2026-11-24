"""proof_controls.py -- controls for the exact machinery.
(1) Mutation test: change one coefficient of G1's B3 by 1/6 -> the all-n identity must FAIL.
(2) B2 (M_2, G1 Theorem B2, HPL form) proved for all n by the same Z-sum identity.
(3) Symbol implementation checks (exact):
    (a) shuffle compatibility  sum_{w in u sha v} S(H_w) == S(H_u) sha S(H_v)  for all |u|+|v| <= 4,
    (b) x -> -x symmetry:  S(H_{-w}) == (-1)^{#nonzero} swap_{1-x <-> 1+x} S(H_w),
    (c) derivative of the last slot: D S(H_{a,w}) == f_a(x) S(H_w) (letters 0,+-1; constant letter 2 drops).
"""
import sys, json, re, pickle, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
sys.setrecursionlimit(10000)
import zsum, hpl
from zsum import K, nn
import b3_allN_proof as P
import mu3_symbolic as M

# ---------------- (1) mutation test
Zmu = pickle.load(open('mu_zsum.pkl', 'rb'))
terms = P.load(P.FILES['1'])
mut = list(terms)
e, w, q = mut[30]
mut[30] = (e, w, q + Fr(1, 6))
ev = P.hpl_side(mut, 0)
muE = {tuple(ww): P.toK(c) for ww, c in Zmu['1']}
diff = zsum.lin((1, ev), (-1, muE))
od = P.hpl_side(mut, 1)
print(f'(1) mutation of term {w} (c^{e}) by +1/6: even identity holds? {len(diff) == 0}; odd identity holds? {len(od) == 0}  (both must be False)')

# ---------------- (2) B2 for all n
G1 = P.G1
F2 = M.mu_formula(2)
Z2 = M.to_zsum(F2)
files2 = {'1': 'hpl_fit_k2_(0,0,0).json', 'z2': 'hpl_fit_k2_(1,0,0).json'}
ok2 = True
for part, fn in files2.items():
    t2 = P.load(fn)
    e2 = P.hpl_side(t2, 0)
    d2 = zsum.lin((1, e2), (-1, Z2[part]))
    o2 = P.hpl_side(t2, 1)
    lo = P.low_order_check(t2, part, lambda nv: M.mu_exact(2, nv))
    print(f'(2) B2 part {part}: {len(t2)} HPL terms; even identity {len(d2) == 0}; odd identity {len(o2) == 0}; low orders {lo}')
    ok2 = ok2 and not d2 and not o2 and lo
print('    B2 (HPL form of M_2) PROVED for all n:', ok2, '; mu_{2,n} has no z3 part:', len(Z2['z3']) == 0)

# ---------------- (3) symbol checks
def shuffle_words(u, v):
    if not u:
        return {tuple(v): 1}
    if not v:
        return {tuple(u): 1}
    out = {}
    for ww, m in shuffle_words(u[1:], v).items():
        out[(u[0],) + ww] = out.get((u[0],) + ww, 0) + m
    for ww, m in shuffle_words(u, v[1:]).items():
        out[(v[0],) + ww] = out.get((v[0],) + ww, 0) + m
    return out

def sym_shuffle(S1, S2):
    out = {}
    for w1, c1 in S1.items():
        for w2, c2 in S2.items():
            for ww, m in shuffle_words(w1, w2).items():
                out[ww] = out.get(ww, 0) + c1*c2*m
    return {k: v for k, v in out.items() if v != 0}

def words(L):
    return [w for w in itertools.product([0, 1, -1], repeat=L) if w[-1] != 0]

ok_a = True; cnt = 0
for L1 in range(1, 4):
    for L2 in range(1, 5 - L1):
        for u in words(L1):
            for v in words(L2):
                lhs = {}
                for ww, m in shuffle_words(u, v).items():
                    for k, c in hpl.symbol(ww).items():
                        lhs[k] = lhs.get(k, 0) + m*c
                lhs = {k: c for k, c in lhs.items() if c != 0}
                if lhs != sym_shuffle(hpl.symbol(u), hpl.symbol(v)):
                    ok_a = False
                cnt += 1
print(f'(3a) shuffle compatibility on {cnt} pairs: {ok_a}')

ok_b = True
swap = {0: 0, 1: 1, 2: 3, 3: 2}
for L in range(1, 5):
    for w in words(L):
        mw = tuple(-t for t in w)
        sg = (-1)**sum(1 for t in w if t != 0)
        S = hpl.symbol(w)
        T = {tuple(swap[i] for i in k): sg*c for k, c in S.items()}
        if hpl.symbol(mw) != T:
            ok_b = False
print('(3b) x -> -x symmetry for all words of weight <= 4:', ok_b)

# (c) last-slot derivative: as a formal combination over the dlog basis {1/x, -1/(1-x), 1/(1+x)}.
# D(S) := sum c (a_1..a_{n-1}) * [dlog a_n]; we compare the coefficient vectors (dlog x, dlog(1-x), dlog(1+x)).
ok_c = True
fa = {0: (1, 0, 0), 1: (0, -1, 0), -1: (0, 0, 1)}   # f_a dx in terms of dlog x, dlog(1-x), dlog(1+x)
for L in range(2, 5):
    for w in words(L):
        S = hpl.symbol(w)
        D = {}
        for k, c in S.items():
            last = k[-1]
            if last == 0:
                continue      # letter 2: constant
            vec = [0, 0, 0]; vec[last - 1] = 1
            for j in range(3):
                if vec[j]:
                    key = (k[:-1], j)
                    D[key] = D.get(key, 0) + c
        D = {k: v for k, v in D.items() if v != 0}
        T = {}
        for k, c in hpl.symbol(w[1:]).items():
            for j, f in enumerate(fa[w[0]]):
                if f:
                    T[(k, j)] = T.get((k, j), 0) + c*f
        T = {k: v for k, v in T.items() if v != 0}
        if D != T:
            ok_c = False
print('(3c) last-slot derivative identity for all words of weight 2..4:', ok_c)
