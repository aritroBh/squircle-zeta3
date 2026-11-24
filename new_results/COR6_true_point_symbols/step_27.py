lib = r'''"""cor6_sym.py -- COR6: TRUE number-level symbols of H_w(c0), constant letters included.

Conventions (pinned by controls, see verify_cor6.py):
  J(a1..an; x): iterated integral, leftmost letter = OUTERMOST:
     J(a1..an;x) = int_0^x dt/(t-a1) J(a2..an; t),  J(empty;x)=1,
     J(0;x) = log x (regularized), J(a;x) = log(1-x/a) for a != 0.
  H_w(x) = (-1)^{#(a_i = +1)} J(w; x)   (HPL, leftmost index = outermost).
  Goncharov symbol recursion (shuffle-regularized):
     S(J(a1..an;x)) = sum_{j=1..n} [l_j] (x) S(J(a1..aj-1, aj+1..an; x)),
     l_j = (a_{j+1}-a_j)/(a_{j-1}-a_j),  a_0 := 0, a_{n+1} := x.
     A zero numerator or denominator factor is DROPPED (regularization);
     if both vanish the term contributes 0.
  Symbols: dict word (tuple of letter vectors) -> Fraction.
"""
from fractions import Fraction as Fr

def ladd(u, v): return tuple(a + b for a, b in zip(u, v))
def lsub(u, v): return tuple(a - b for a, b in zip(u, v))
def lneg(u): return tuple(-a for a in u)

def sadd(*S):
    out = {}
    for s in S:
        for w, c in s.items():
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def sscale(c, s):
    return {w: v * c for w, v in s.items() if v * c}

def shuffles(u, v):
    if not u: return [v]
    if not v: return [u]
    return ([(u[0],) + w for w in shuffles(u[1:], v)]
            + [(v[0],) + w for w in shuffles(u, v[1:])])

def sshuffle(s, t):
    out = {}
    for u, cu in s.items():
        for v, cv in t.items():
            for w in shuffles(u, v):
                out[w] = out.get(w, 0) + cu * cv
    return {w: c for w, c in out.items() if c}

# ---------------- Goncharov symbol ----------------
def sym_J(word, x, vec, a0=0, prepend=True, inv=False):
    """Symbol of J(word; x). vec(el) -> letter vector (tuple of Fr); el != 0 always.
    Elements need == and -. Memoized per (vec id) by caller wrapping."""
    memo = {}
    def rec(word):
        n = len(word)
        if n == 0: return {}
        if word in memo: return memo[word]
        res = {}
        for j in range(n):
            p = a0 if j == 0 else word[j - 1]
            q = x if j == n - 1 else word[j + 1]
            wj = word[j]
            num = None if q == wj else q - wj
            den = None if p == wj else p - wj
            if num is None and den is None:
                continue
            vn = None if num is None else vec(num)
            vd = None if den is None else vec(den)
            if inv: vn, vd = vd, vn
            if vn is None and vd is None:
                continue
            lv = lneg(vd) if vn is None else (vn if vd is None else lsub(vn, vd))
            for w2, c2 in rec(word[:j] + word[j + 1:]).items():
                w3 = (lv,) + w2 if prepend else w2 + (lv,)
                res[w3] = res.get(w3, 0) + c2
        res = {w: c for w, c in res.items() if c}
        memo[word] = res
        return res
    return rec(tuple(word))

def sym_H(word, x, vec, **kw):
    """Symbol of HPL H_word(x)."""
    s = -1 if sum(1 for a in word if a == 1) % 2 else 1
    return sscale(s, sym_J(tuple(word), x, vec, **kw))

# ---------------- Lie / quotient machinery ----------------
def expand_word(w, dim):
    terms = [((), Fr(1))]
    for l in w:
        terms = [(p + (ej,), c * l[j]) for p, c in terms for j in range(dim)
                 if l[j] for ej in [tuple(Fr(1) if i == j else Fr(0) for i in range(dim))]]
    out = {}
    for p, c in terms:
        out[p] = out.get(p, 0) + c
    return {p: c for p, c in out.items() if c}

def expand_sym(S, dim):
    out = {}
    for w, cf in S.items():
        for p, c in expand_word(w, dim).items():
            out[p] = out.get(p, 0) + cf * c
    return {p: c for p, c in out.items() if c}

def rref(rows):
    rows = [list(r) for r in rows]; R = 0; piv = []
    for j in range(len(rows[0])):
        p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
        if p is None: continue
        rows[R], rows[p] = rows[p], rows[R]
        q = rows[R][j]; rows[R] = [x / q for x in rows[R]]
        for i in range(len(rows)):
            if i != R and rows[i][j]:
                rows[i] = [a - b * rows[i][j] for a, b in zip(rows[i], rows[R])]
        piv.append(j); R += 1
        if R == len(rows): break
    return rows[:R], piv

class W4Quotient:
    """Q^4 = T^4(V)/D^4 for V = Q^dim: weight-4 symbols modulo shuffle products."""
    def __init__(self, dim):
        from itertools import product as ip
        self.dim = dim
        self.basis_letters = [tuple(Fr(1) if i == j else Fr(0) for i in range(dim)) for j in range(dim)]
        self.W4 = list(ip(self.basis_letters, repeat=4))
        W1 = self.basis_letters
        W2 = list(ip(self.basis_letters, repeat=2))
        W3 = list(ip(self.basis_letters, repeat=3))
        D = []
        for u in W1:
            for v in W3: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), dim))
        for u in W2:
            for v in W2: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), dim))
        def vec_of(S): return tuple(S.get(w, Fr(0)) for w in self.W4)
        self.vec_of = vec_of
        self.Dmat, self.piv = rref([vec_of(s) for s in D])
        self.free = [j for j in range(len(self.W4)) if j not in self.piv]
        self.free_words = [self.W4[j] for j in self.free]
    def coords(self, S):
        v = list(self.vec_of(expand_sym(S, self.dim)))
        for i, pj in enumerate(self.piv):
            if v[pj]:
                q = v[pj]
                v = [a - q * b for a, b in zip(v, self.Dmat[i])]
        assert all(v[j] == 0 for j in range(len(self.W4)) if j not in self.free), "not integrable/decomposable?"
        return tuple(v[j] for j in self.free)

# ---------------- delta_{2,2} (COR4 convention: (l1^l2)^(l3^l4), antisymmetrised) ----------------
def wedge(u, v, dim):
    out = {}
    for i in range(dim):
        for j in range(i + 1, dim):
            c = u[i] * v[j] - u[j] * v[i]
            if c: out[(i, j)] = out.get((i, j), 0) + c
    return out

def delta22(S, dim):
    """(a,b,c,d) -> (a^b)^(c^d), antisymmetrised within and between pairs.
    Kills Li4 and all shuffle products. Returns dict ((ij),(kl)) -> Fr, (ij)<=(kl)."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2: continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                out[key] = out.get(key, 0) + cf * c1 * c2 * sgn
    return {k: v for k, v in out.items() if v}

# ---------------- planner's rho (Lie coordinates) ----------------
def rho(T):
    """T: dict basis-word (tuple of basis-letter indices) -> Fr.
    Returns dict word -> Fr via the planner's recursion."""
    memo = {}
    def r(word):
        if word in memo: return memo[word]
        n = len(word)
        if n == 1: out = {word: Fr(1)}
        else:
            out = {}
            for k, v in r(word[:-1]).items():
                out[k + (word[-1],)] = out.get(k + (word[-1],), 0) + Fr(n - 1, n) * v
            for k, v in r(word[1:]).items():
                out[k + (word[0],)] = out.get(k + (word[0],), 0) - Fr(n - 1, n) * v
        memo[word] = out
        return out
    tot = {}
    for w, co in T.items():
        for k, v in r(w).items():
            tot[k] = tot.get(k, 0) + co * v
    return {k: v for k, v in tot.items() if v != 0}
'''
with open("cor6_sym.py", "w") as f:
    f.write(lib)
print("ok")
