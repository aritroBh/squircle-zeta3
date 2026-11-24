"""g8_M3_delta22.py -- Proposition 5 (sec_g8_M3alpha.tex): the (2,2) cobracket of Omega(c), and the two-term identity.
Exact arithmetic (python Fractions); written from the definitions of the section.  Reads the 56 weight-four terms of R_3
from g8_terms_k3_part1.json (the data of the Table in the section; the weight of a term is the length of its word).

Conventions (section g8:sec:M3alpha).
  V = Q-span of the letters C=[c], A=[1-c], B=[1+c].   S(H_{a1..an}(c)) = g(an) (x) ... (x) g(a1), g(0)=C, g(1)=-A, g(-1)=B.
  A coefficient q c^e is a scalar (Laurent polynomial in c).
  delta22 (x1 (x) x2 (x) x3 (x) x4) = (x1 ^ x2) ^ (x3 ^ x4) in Lambda^2(Lambda^2 V).
Part 1: S(Omega) and delta22(Omega).   Part 2: controls.   Part 3: one term cannot.   Part 4: two terms can (mod [2]), and the
exact residual when the constant letter [2] is kept.
"""
import json, re, itertools, random
from fractions import Fraction as Fr
from collections import defaultdict

# ------------------------------------------------ generic linear algebra on tensors with Laurent coefficients
# a "Laurent polynomial" is a dict {e: Fraction}; a vector is {key: Laurent}
def ladd(p, q, s=1):
    out = dict(p)
    for e, c in q.items():
        v = out.get(e, 0) + s*c
        if v == 0: out.pop(e, None)
        else: out[e] = v
    return out
def lscale(p, c):
    return {e: c*v for e, v in p.items() if c*v != 0}
def lshift(p, k):
    return {e + k: v for e, v in p.items()}
def vadd(V, key, lau):
    if not lau: return
    cur = V.get(key, {})
    new = ladd(cur, lau)
    if new: V[key] = new
    else: V.pop(key, None)

def lstr(p):
    if not p: return '0'
    return ' + '.join(f"({c})c^{e}" for e, c in sorted(p.items()))

# ------------------------------------------------ part 1: symbol of H_w
G = {0: {'C': Fr(1)}, 1: {'A': Fr(-1)}, -1: {'B': Fr(1)}}      # g(a) as a vector in V
def tensor_expand(letters_vectors):
    """product of vectors -> dict {tuple of basis letters: coefficient}"""
    out = {(): Fr(1)}
    for vec in letters_vectors:
        new = {}
        for t, c in out.items():
            for l, d in vec.items():
                new[t + (l,)] = new.get(t + (l,), 0) + c*d
        out = new
    return out

def symbol_H(w):
    """S(H_w(c)) as {tuple: Fraction}; tensor g(a_n) ... g(a_1)."""
    return tensor_expand([G[a] for a in reversed(w)])

terms = []
for nm, v in json.load(open('g8_terms_k3_part1.json')):
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
    terms.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
omega = [t for t in terms if len(t[1]) == 4]
print(f"R_3 has {len(terms)} terms; {len(omega)} of weight four (Omega); weights 2,3,4: "
      f"{sum(1 for t in terms if len(t[1])==2)}, {sum(1 for t in terms if len(t[1])==3)}, {len(omega)}")

def symbol_of_terms(ts):
    S = {}
    for e, w, q in ts:
        for t, c in symbol_H(w).items():
            vadd(S, t, {e: q*c})
    return S

S_om = symbol_of_terms(omega)
print(f"S(Omega): {len(S_om)} pure tensors with nonzero Laurent coefficient (letters C=[c], A=[1-c], B=[1+c]); exponents of c:",
      sorted({e for v in S_om.values() for e in v}))

ORDER = {'C': 0, 'A': 1, 'B': 2, 'T': 3, 'X': 4, 'U': 5, 'Y': 6, 'V': 7, 'W': 8}   # any fixed total order on basis letters

def wedge(a, b):
    """a ^ b in Lambda^2 V for basis letters -> (sign, (i<j)) or None"""
    if a == b: return None
    return (1, (a, b)) if ORDER[a] < ORDER[b] else (-1, (b, a))

def delta22(S):
    """sum over tensors: (x1^x2)^(x3^x4); returns {(pair1, pair2): Laurent} with pair1 < pair2 canonical."""
    out = {}
    for t, lau in S.items():
        w1 = wedge(t[0], t[1]); w2 = wedge(t[2], t[3])
        if w1 is None or w2 is None: continue
        s = w1[0]*w2[0]
        p1, p2 = w1[1], w2[1]
        key_order = lambda p: (ORDER[p[0]], ORDER[p[1]])
        if p1 == p2: continue
        if key_order(p1) > key_order(p2):
            p1, p2 = p2, p1; s = -s
        vadd(out, (p1, p2), lscale(lau, s))
    return out

def show22(D):
    for (p1, p2), lau in sorted(D.items()):
        print(f"    ([{p1[0]}]^[{p1[1]}]) ^ ([{p2[0]}]^[{p2[1]}])  with coefficient  {lstr(lau)}")

D_om = delta22(S_om)
print("delta22(Omega) =")
show22(D_om)
# the statement of Proposition 5: (1 - 1/c) (C^A)^(A^B) - (1 + 1/c) (C^B)^(A^B)
target = {(( 'C', 'A'), ('A', 'B')): {0: Fr(1), -1: Fr(-1)}, (('C', 'B'), ('A', 'B')): {0: Fr(-1), -1: Fr(-1)}}
print("equals the displayed formula of Proposition 5:", D_om == target)
print("it is not zero (two nonzero components, coefficients 1 - 1/c and -(1 + 1/c) are nonzero Laurent polynomials)")

# ------------------------------------------------ part 2: controls (delta22 kills products and Li_4)
def shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + s for s in shuffle(u[1:], v)] + [(v[0],) + s for s in shuffle(u, v[1:])]
random.seed(1)
letters = ['C', 'A', 'B']
bad = 0; cnt = 0
for _ in range(300):
    for (p, q) in ((1, 3), (2, 2), (3, 1)):
        u = tuple(random.choice(letters) for _ in range(p)); v = tuple(random.choice(letters) for _ in range(q))
        S = {}
        for s in shuffle(u, v): vadd(S, s, {0: Fr(1)})
        cnt += 1
        if delta22(S): bad += 1
print(f"control 1: delta22 of {cnt} random shuffle products (1,3),(2,2),(3,1): nonzero in {bad}")
# exhaustive: all products of letters
bad = 0; cnt = 0
for p in (1, 2, 3):
    for u in itertools.product(letters, repeat=p):
        for v in itertools.product(letters, repeat=4 - p):
            S = {}
            for s in shuffle(u, v): vadd(S, s, {0: Fr(1)})
            cnt += 1
            if delta22(S): bad += 1
print(f"control 2: delta22 of ALL {cnt} shuffle products of words over {{C,A,B}} of total length 4: nonzero in {bad}")
# Li_4(r): -[1-r](x)[r]^3 with r, 1-r any two letters (here all ordered pairs of basis letters, also repeated)
bad = 0
for a, b in itertools.product(['C', 'A', 'B', 'T'], repeat=2):
    S = {(a, b, b, b): {0: Fr(-1)}}
    if delta22(S): bad += 1
print("control 3: delta22 of -[a](x)[b](x)[b](x)[b] (symbol of Li_4(r) with [a]=[1-r], [b]=[r]) nonzero in", bad, "of 16 pairs")
# positive control: a symbol that is not classical
S = {('C', 'A', 'A', 'B'): {0: Fr(1)}}
print("control 4 (positive): delta22 of C(x)A(x)A(x)B is nonzero:", bool(delta22(S)))

# ------------------------------------------------ part 3: one term cannot
A_ = D_om[(('C', 'A'), ('A', 'B'))]; B_ = D_om[(('C', 'B'), ('A', 'B'))]
print("part 3: the two components are", lstr(A_), "and", lstr(B_))
def lmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            out[e1 + e2] = out.get(e1 + e2, 0) + c1*c2
    return {e: v for e, v in out.items() if v != 0}
lhs = lmul(A_, {0: Fr(1), 1: Fr(1)}); rhs = lmul(B_, {0: Fr(1), 1: Fr(-1)})
print("   A (1+c) =", lstr(lhs), "  and  B (1-c) =", lstr(rhs), "; equal:", lhs == rhs, " (so A/B = (1-c)/(1+c), not a constant)")

# ------------------------------------------------ part 4: two-variable symbols of Li_{a,b}(x,y)
# basis letters X=[x], U=[1-x], Y=[y], V=[1-y], W=[1-xy]; Z=[xy]=X+Y.   Vectors in this 5-dim space (constant coefficients).
def vec(**kw): return {k: Fr(v) for k, v in kw.items() if v}
Zv = vec(X=1, Y=1)
def tens(*vecs):          # tensor product of vectors -> {tuple: Fraction}
    return tensor_expand(list(vecs))
def tsum(*Ts):
    out = {}
    for sgn, T in Ts:
        for k, v in T.items():
            out[k] = out.get(k, 0) + sgn*v
    return {k: v for k, v in out.items() if v != 0}
def tapp(T, vecw):        # T (x) vecw
    out = {}
    for k, v in T.items():
        for l, d in vecw.items():
            out[k + (l,)] = out.get(k + (l,), 0) + v*d
    return {k: v for k, v in out.items() if v != 0}
Wv, Uv, Vv, Xv, Yv = vec(W=1), vec(U=1), vec(V=1), vec(X=1), vec(Y=1)
# S(Li_1(z)) = -[1-z];  S(Li_2(z)) = -[1-z](x)[z]
S11 = tsum((1, tapp({('W',): Fr(1)}, vec(Y=1, U=1, V=-1))), (1, tapp({('U',): Fr(1)}, Vv)))
S12 = tsum((1, tens(Wv, Zv, Uv)), (1, tapp(S11, Yv)))
S21 = tsum((1, tapp(S11, Xv)), (1, tapp(tens(Wv, Zv), vec(Y=1, V=-1))), (1, tens(Uv, Xv, Vv)))
S22 = tsum((1, tapp(S12, Xv)), (1, tapp(S21, Yv)))
S31 = tsum((1, tapp(S21, Xv)), (1, tapp(tens(Wv, Zv, Zv), vec(Y=1, V=-1))), (1, tens(Uv, Xv, Xv, Vv)))
# (dLi_{3,1} = Li_{2,1} dlog x - Li_3(xy) dlog(y/(y-1)) - Li_3(x) dlog(1-y))
S31 = tsum((1, tapp(S21, Xv)), (1, tapp(tens(Wv, Zv, Zv), vec(Y=1, V=-1))), (1, tens(Uv, Xv, Xv, Vv)))

def expand_Z(T):          # Z letter already expanded in Zv; nothing to do (kept for clarity)
    return T

print("two-variable symbols: Li_{1,1}:", len(S11), "terms, Li_{2,2}:", len(S22), "tensors, Li_{3,1}:", len(S31), "tensors (letters X,U,Y,V,W, Z=X+Y expanded)")

def d22_const(T):
    return delta22({k: {0: v} for k, v in T.items()})

def specialise(T, mp_, drop_last_const=False, const_letters=()):
    """apply a letter map  mp_: letter -> vector in letters C,A,B,T ; returns {tuple: Fraction}."""
    out = {}
    for k, v in T.items():
        for t, c in tensor_expand([mp_[l] for l in k]).items():
            out[t] = out.get(t, 0) + v*c
    out = {k: v for k, v in out.items() if v != 0}
    if drop_last_const:
        out = {k: v for k, v in out.items() if k[-1] not in const_letters}
    return out

# x = 1/2:  X = U = [1/2] = -T   (T = [2] constant letter).  y = 1-c and y = 1+c.
def mapping(case, keepT):
    t = -1 if keepT else 0
    Xm = {'T': Fr(t)} if keepT else {}
    if case == 'minus':     # y = 1 - c : Y=[1-c]=A, V=[c]=C, xy=(1-c)/2: Z = A - T, 1-xy=(1+c)/2: W = B - T
        Y_ = {'A': Fr(1)}; V_ = {'C': Fr(1)}; W_ = {'B': Fr(1)}
        if keepT: W_['T'] = Fr(-1)
        Z_ = {'A': Fr(1)}
        if keepT: Z_['T'] = Fr(-1)
    else:                   # y = 1 + c : Y=B, V=[-c]=C, xy=(1+c)/2: Z = B - T, 1-xy=(1-c)/2: W = A - T
        Y_ = {'B': Fr(1)}; V_ = {'C': Fr(1)}; W_ = {'A': Fr(1)}
        if keepT: W_['T'] = Fr(-1)
        Z_ = {'B': Fr(1)}
        if keepT: Z_['T'] = Fr(-1)
    return {'X': dict(Xm), 'U': dict(Xm), 'Y': Y_, 'V': V_, 'W': W_}

def Li22_symbol(case, keepT):
    mp_ = mapping(case, keepT)
    return specialise(S22, mp_)

# generic delta22 of Li_{2,2}, Li_{3,1} (formulas of COR8BR): print
D22g = d22_const(S22)
print("delta22(Li_{2,2}(x,y)) generic (5 letters X,U,Y,V,W; Z=X+Y):")
show22(D22g)
D31g = d22_const(S31)
print("delta22(Li_{3,1}(x,y)) generic:")
show22(D31g)

for keepT in (False, True):
    tag = "constant letter [2] kept as a formal letter, no letter dropped (NOT used in the section)" if keepT else "iterated-differential symbol: a constant has no letter (X, U -> 0)"
    Sm = Li22_symbol('minus', keepT); Sp = Li22_symbol('plus', keepT)
    Dm = delta22({k: {0: v} for k, v in Sm.items()})
    Dp = delta22({k: {0: v} for k, v in Sp.items()})
    # DT = ((1-c)/c) Li22(1/2,1-c) - ((1+c)/c) Li22(1/2,1+c) ;  (1-c)/c = c^-1 - 1 ;  (1+c)/c = c^-1 + 1
    DT = {}
    for key, lau in Dm.items(): vadd(DT, key, lmul(lau, {-1: Fr(1), 0: Fr(-1)}))
    for key, lau in Dp.items(): vadd(DT, key, lmul(lau, {-1: Fr(-1), 0: Fr(-1)}))
    print(f"\n[{tag}]  delta22 of Li_{{2,2}}(1/2,1-c):")
    show22(Dm)
    print(f"[{tag}]  delta22 of Li_{{2,2}}(1/2,1+c):")
    show22(Dp)
    print(f"[{tag}]  delta22(DT), DT = ((1-c)/c) Li_{{2,2}}(1/2,1-c) - ((1+c)/c) Li_{{2,2}}(1/2,1+c):")
    show22(DT)
    # difference with delta22(Omega)
    R = {}
    for key, lau in DT.items(): vadd(R, key, lau)
    for key, lau in D_om.items(): vadd(R, key, lscale(lau, Fr(-1)))
    print(f"[{tag}]  delta22(DT) - delta22(Omega) has {len(R)} nonzero components" + (":" if R else " (zero)"))
    show22(R)
