"""g8b_check.py -- checks for sec_g8_M3func.tex (Proposition 5, function level).

Conventions of the section.  Alphabet Sigma (finite set of points), base point z0.
  Hyperlogarithm  L_{s1 s2 ... sn}(z) = int_{z0}^z dt/(t-s1) L_{s2...sn}(t),  L_empty = 1.
  rho(s1...sn) = e_{sn} (x) ... (x) e_{s1}    (outermost integration letter s1 LAST), e_s = [z - s] mod constants.
  For F in the C-span of the L_w:  S(F) = sum_w c_w rho(w), graded by length; S_k = degree-k part.
  delta22(x1 x2 x3 x4) = (x1^x2)^(x3^x4).

Part 1 (exact).  T_4(M_3) from the data g8_terms_k3_part1.json / g8_terms_k3_partz2.json: word lengths; the top part of
        H_w is s_w L_w with s_w = (-1)^{#letters 1}; rho of the 56 weight-four terms against the 24-row table of
        sec_g8_M3alpha.tex (copied by hand below), and delta22 against the formula of Proposition 5.
Part 2 (exact, sympy).  S_4 of branches of Li_n(g) (n = 3..6), of products (Li_2 Li_2, Li_3 log, Li_2 Li_3, log^4, Li_1^2 Li_2)
        for sample g with letters in {0,1,-1}, with the base-point constants P_k(z0) kept as formal symbols; delta22 = 0.
Part 3 (numerical, mpmath, 50 digits).  (a) Li_4(g(z)) and Li_5(g(z)), g = (1-z)/(1+z), equal their expansion
        sum_w c_w L_w(z) with base point z0 = 1/2 (Lemma g8b:lem:classical);  (b) the base-point formula
        H_w(z) = sum_{w=uv} s_u L_u(z) H_v(z0) for words of the table (Lemma g8b:lem:H).
"""
import json, re, itertools, sys
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

HERE = __file__.rsplit('/', 1)[0] if '/' in __file__ else '.'

# ------------------------------------------------------------------ tensors and delta22 (generic coefficients)
def tadd(T, key, val):
    v = T.get(key, 0) + val
    if v == 0:
        T.pop(key, None)
    else:
        T[key] = v

def wedge(a, b, order):
    if a == b:
        return None
    return (1, (a, b)) if order[a] < order[b] else (-1, (b, a))

def delta22(T, order, simplify=lambda v: v):
    out = {}
    for t, val in T.items():
        w1 = wedge(t[0], t[1], order); w2 = wedge(t[2], t[3], order)
        if w1 is None or w2 is None:
            continue
        s = w1[0]*w2[0]; p1, p2 = w1[1], w2[1]
        if p1 == p2:
            continue
        k1 = (order[p1[0]], order[p1[1]]); k2 = (order[p2[0]], order[p2[1]])
        if k1 > k2:
            p1, p2, s = p2, p1, -s
        tadd(out, (p1, p2), s*val)
    return {k: simplify(v) for k, v in out.items() if simplify(v) != 0}

def shuffle_words(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + s for s in shuffle_words(u[1:], v)] + [(v[0],) + s for s in shuffle_words(u, v[1:])]

def graded_shuffle(A, B):
    """A, B: {tuple: coeff} (all degrees) -> shuffle product."""
    out = {}
    for a, x in A.items():
        for b, y in B.items():
            for s in shuffle_words(a, b):
                tadd(out, s, x*y)
    return out

# ------------------------------------------------------------------ Part 1
print("=== Part 1: T_4(M_3) from the data ===")
def read_terms(fn):
    out = []
    for nm, v in json.load(open(f"{HERE}/{fn}")):
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        if m is None:
            out.append((0, (), Fr(v)))          # the constant term "1"
        else:
            out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out
R3 = read_terms('g8_terms_k3_part1.json')
Psi = read_terms('g8_terms_k3_partz2.json')
from collections import Counter
print("R_3: number of terms by word length:", dict(sorted(Counter(len(w) for _, w, _ in R3).items())))
print("Psi: number of terms by word length:", dict(sorted(Counter(len(w) for _, w, _ in Psi).items())),
      "(length 0 = the constant 1/2)")
print("Phi = 1/2 - (1-c^2)/(4c) (H_1 + H_{-1}): length <= 1")
print("all words regular (last letter nonzero):", all(w[-1] != 0 for _, w, _ in R3 + Psi if w))
# top part:  H_w = s_w L_w + shorter;  rho(w) = e_{a_n} ... e_{a_1}
NAME = {0: 'x', 1: 'y', -1: 'z'}          # e_0 = [c], e_1 = [c-1] = [1-c], e_{-1} = [c+1]   (mod constants)
T4 = {}
for e, w, q in R3:
    if len(w) != 4:
        continue
    s = (-1)**sum(1 for a in w if a == 1)
    key = tuple(NAME[a] for a in reversed(w))
    cur = dict(T4.get(key, {}))
    cur[e] = cur.get(e, 0) + s*q
    cur = {k: v for k, v in cur.items() if v != 0}
    if cur: T4[key] = cur
    else: T4.pop(key, None)
print("rho(T_4(M_3)) has", len(T4), "pure tensors")
# the table of sec_g8_M3alpha.tex (coefficient as {exponent: value}), copied by hand
h = Fr(1, 2)
TABLE = {
 'yxxx': {0: 1, 1: -1}, 'yxxy': {-1: -h, 1: h}, 'yxxz': {-1: h, 0: 1, 1: h}, 'yxyx': {-1: -h, 1: h},
 'yxyz': {-1: h, 1: -h}, 'yxzx': {-1: h, 0: 1, 1: h}, 'yxzy': {-1: -h, 0: 1, 1: -h}, 'yzxx': {0: 1, 1: 1},
 'yzxy': {-1: -h, 0: 1, 1: -h}, 'yzxz': {-1: h, 1: -h}, 'yzyx': {0: 1, 1: -1}, 'yzyz': {0: 1, 1: 1},
 'zxxx': {0: 1, 1: 1}, 'zxxy': {-1: -h, 0: 1, 1: -h}, 'zxxz': {-1: h, 1: -h}, 'zxyx': {-1: -h, 0: 1, 1: -h},
 'zxyz': {-1: h, 0: 1, 1: h}, 'zxzx': {-1: h, 1: -h}, 'zxzy': {-1: -h, 1: h}, 'zyxx': {0: 1, 1: -1},
 'zyxy': {-1: -h, 1: h}, 'zyxz': {-1: h, 0: 1, 1: h}, 'zyzx': {0: 1, 1: 1}, 'zyzy': {0: 1, 1: -1}}
TABLE = {tuple(k): {e: Fr(v) for e, v in d.items()} for k, d in TABLE.items()}
print("rho(T_4(M_3)) equals the 24-row table of sec_g8_M3alpha.tex:", T4 == TABLE)
c = sp.Symbol('c')
T4s = {k: sum(sp.Rational(v.numerator, v.denominator)*c**e for e, v in d.items()) for k, d in T4.items()}
ORD = {'x': 0, 'y': 1, 'z': 2}
D = delta22(T4s, ORD, simplify=lambda v: sp.simplify(v))
for k, v in sorted(D.items()):
    print("   ", k, ":", sp.factor(v))
target = {(('x', 'y'), ('y', 'z')): 1 - 1/c, (('x', 'z'), ('y', 'z')): -(1 + 1/c)}
print("delta22(T_4(M_3)) equals (1-1/c) p^r - (1+1/c) q^r:",
      set(D) == set(target) and all(sp.simplify(D[k] - target[k]) == 0 for k in target))

# ------------------------------------------------------------------ Part 2
print("\n=== Part 2: S_4 of classical functions (formal base-point constants) ===")
ORD2 = {0: 0, 1: 1, -1: 2}
def cls(factors):
    """[g] for g = const * prod (z - s)^m  given as {s: m}: vector {s: m}."""
    return {s: m for s, m in factors.items() if m}
def vec_to_T(v):
    return {(s,): sp.Integer(m) for s, m in v.items()}
def tensor_append(T, v):
    out = {}
    for k, a in T.items():
        for s, m in v.items():
            tadd(out, k + (s,), a*m)
    return out
def S_log(name, g):
    """full graded S of a branch of log g: constant + [g]."""
    out = {(): sp.Symbol(f'{name}0')}
    for k, a in vec_to_T(g).items():
        tadd(out, k, a)
    return out
def S_Li(name, n, g, one_minus_g):
    """full graded S of the n-th member of a chain P_1..P_n with exp(-P_1)=1-g, P_k' = P_{k-1} g'/g.
    Base-point constants P_k(z0) are formal symbols name_k."""
    S = {(): sp.Symbol(f'{name}_1')}
    for k, a in vec_to_T(one_minus_g).items():
        tadd(S, k, -a)
    for kk in range(2, n + 1):
        S2 = tensor_append(S, g)
        tadd(S2, (), sp.Symbol(f'{name}_{kk}'))
        S = S2
    return S
def deg(T, k):
    return {t: v for t, v in T.items() if len(t) == k}
samples = {   # name: ([g], [1-g])   over the letters 0, 1, -1
    '(1-z)/(1+z)': (cls({1: 1, -1: -1}), cls({0: 1, -1: -1})),
    'z^2':         (cls({0: 2}), cls({1: 1, -1: 1})),
    '4z/(1+z)^2':  (cls({0: 1, -1: -2}), cls({1: 2, -1: -2})),
    '-z':          (cls({0: 1}), cls({-1: 1})),
    'z':           (cls({0: 1}), cls({1: 1})),
}
ok = True
for nm, (g, omg) in samples.items():
    for n in range(3, 7):
        S = S_Li('P', n, g, omg)
        S4 = deg(S, 4)
        d = delta22(S4, ORD2, simplify=sp.expand)
        if n == 4:
            # compare with -[1-g] (x) [g]^3
            ref = {}
            for t, a in vec_to_T(omg).items():
                ref[t] = -a
            for _ in range(3): ref = tensor_append(ref, g)
            same = all(sp.expand(S4.get(k, 0) - ref.get(k, 0)) == 0 for k in set(S4) | set(ref))
        elif n > 4:
            ref = {(): sp.Symbol(f'P_{n-4}')}
            for _ in range(4): ref = tensor_append(ref, g)
            same = all(sp.expand(S4.get(k, 0) - ref.get(k, 0)) == 0 for k in set(S4) | set(ref))
        else:
            same = (S4 == {})
        ok = ok and same and not d
        print(f"  Li_{n}(g), g = {nm}: S_4 has {len(S4)} tensors; equals the formula of Lemma g8b:lem:classical: {same};"
              f" delta22(S_4) = 0: {not d}")
gA, oA = samples['(1-z)/(1+z)']; gB, oB = samples['z^2']; gC, oC = samples['4z/(1+z)^2']
prods = {
    'Li_2(gA) Li_2(gB)': [S_Li('A', 2, gA, oA), S_Li('B', 2, gB, oB)],
    'Li_3(gA) log(gC)': [S_Li('A', 3, gA, oA), S_log('L', gC)],
    'Li_2(gC) Li_3(gB)  (length 5)': [S_Li('C', 2, gC, oC), S_Li('B', 3, gB, oB)],
    'log(gA)^4': [S_log('L', gA)]*4,
    'Li_1(gA) Li_1(gB) Li_2(gC)': [S_Li('A', 1, gA, oA), S_Li('B', 1, gB, oB), S_Li('C', 2, gC, oC)],
    'Li_4(gB) Li_2(gA)  (length 6)': [S_Li('B', 4, gB, oB), S_Li('A', 2, gA, oA)],
}
for nm, fs in prods.items():
    S = fs[0]
    for f in fs[1:]:
        S = graded_shuffle(S, f)
    S4 = deg(S, 4)
    d = delta22(S4, ORD2, simplify=sp.expand)
    ok = ok and not d
    print(f"  product {nm}: S_4 has {len(S4)} tensors; delta22(S_4) = 0: {not d}")
print("Part 2 all as claimed:", ok)

# ------------------------------------------------------------------ Part 3 (numerical)
print("\n=== Part 3: numerical checks of the expansions (mpmath, 50 digits) ===")
mp.mp.dps = 50
NT = 160                      # Taylor order around z0; |z - z0| <= 0.2, radius 0.5: error ~ 0.4^160
z0 = mp.mpf(1)/2
def series_inv(s):            # Taylor coefficients of 1/(z - s) around z0
    a = z0 - s
    return [(-1)**k/a**(k + 1) for k in range(NT)]
INV = {s: series_inv(s) for s in (0, 1, -1)}
def mul(A, B):
    out = [mp.mpf(0)]*NT
    for i, a in enumerate(A):
        if a == 0: continue
        for j in range(NT - i):
            out[i + j] += a*B[j]
    return out
def integ(A):                 # antiderivative vanishing at z0
    return [mp.mpf(0)] + [A[k]/(k + 1) for k in range(NT - 1)]
Lcache = {(): [mp.mpf(1)] + [mp.mpf(0)]*(NT - 1)}
def Lser(w):
    if w not in Lcache:
        Lcache[w] = integ(mul(INV[w[0]], Lser(w[1:])))
    return Lcache[w]
def evalser(A, z):
    t = z - z0
    return mp.fsum(a*t**k for k, a in enumerate(A))
def words_of(T):              # tensor e_{a_n}...e_{a_1} -> word a_1...a_n
    return {tuple(reversed(k)): v for k, v in T.items()}
g = lambda z: (1 - z)/(1 + z)
gv, omgv = samples['(1-z)/(1+z)']
worst = mp.mpf(0)
for n in (4, 5):
    S = S_Li('P', n, gv, omgv)
    vals = {f'P_{k}': mp.polylog(k, g(z0)) for k in range(1, n + 1)}
    for z in (mp.mpf('0.31'), mp.mpf('0.62'), mp.mpf('0.5') + mp.mpf('0.15')*1j):
        tot = mp.mpf(0)
        for w, coef in words_of(S).items():
            expr = sp.expand(coef)
            # coefficient: polynomial in the formal symbols P_k (degree <= 1)
            cnum = mp.mpf(str(expr.subs({sp.Symbol(k): 0 for k in vals}))) if expr.free_symbols else mp.mpf(str(expr))
            for sym in expr.free_symbols:
                cnum += mp.mpf(str(expr.coeff(sym))) * vals[str(sym)]
            tot += cnum*evalser(Lser(w), z)
        diff = abs(tot - mp.polylog(n, g(z)))
        worst = max(worst, diff)
        print(f"  Li_{n}((1-z)/(1+z)) at z = {mp.nstr(z, 4)}: |expansion - Li_{n}| = {mp.nstr(diff, 3)}")
# (b) H_w via its series at 0, against the base-point formula
NH = 400
def H_series(w):              # Taylor coefficients at 0 of H_w (regular word), exact rationals
    if not w:
        return [Fr(1)] + [Fr(0)]*(NH - 1)
    inner = H_series(w[1:])
    a = w[0]
    if a == 0:
        f_times = inner[1:] + [Fr(0)]            # divide by t (inner has no constant term)
    else:
        sgn = 1 if a == 1 else -1                # 1/(1-t) = sum t^k ; 1/(1+t) = sum (-t)^k
        f_times = [Fr(0)]*NH
        acc = Fr(0)
        for k in range(NH):
            acc = acc*sgn + inner[k]
            f_times[k] = acc
    return [Fr(0)] + [f_times[k]/(k + 1) for k in range(NH - 1)]
Hc = {}
def Hval(w, z):
    if w not in Hc:
        Hc[w] = H_series(w)
    return mp.fsum(mp.mpf(a.numerator)/a.denominator*z**k for k, a in enumerate(Hc[w]) if a)
testwords = [(0, 1, 0, -1), (-1, 1, -1, 1), (1, 0, -1, 1), (0, 0, 0, 1), (0, -1, 1, -1)]
for w in testwords:
    for z in (mp.mpf('0.31'), mp.mpf('0.62')):
        tot = mp.mpf(0)
        for k in range(len(w) + 1):
            u, v = w[:k], w[k:]
            su = (-1)**sum(1 for a in u if a == 1)
            tot += su*evalser(Lser(u), z)*(Hval(v, z0) if v else 1)
        diff = abs(tot - Hval(w, z))
        worst = max(worst, diff)
        print(f"  H_{w} at z = {mp.nstr(z, 3)}: |sum_(w=uv) s_u L_u(z) H_v(1/2) - H_w(z)| = {mp.nstr(diff, 3)}")
print("largest difference in Part 3:", mp.nstr(worst, 3), "; below 1e-40:", worst < mp.mpf('1e-40'))
