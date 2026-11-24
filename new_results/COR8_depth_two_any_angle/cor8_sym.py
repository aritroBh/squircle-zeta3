"""cor8_sym.py -- exact symbol algebra for COR8.

Letters: exponent vectors over basis {c, 1-c, 1+c} (indices 0,1,2); constants killed
(function-field letters k(c)^*/k^*).  Extra bookkeeping letters may be appended
(indices >= 3) to certify independence of physical letters (the epsilon trick).

Coefficients: Laurent polynomials in c: dict e -> Fraction (e in Z).
A symbol is a dict word(tuple of letter-vectors) -> Laurent poly.

TRUE symbols (constant coefficients) of two-variable MPLs, derived from the
differential recursion with 1-forms converted to exact dlog form
  x/(1-x) dlog x = -dlog(1-x),   dh/(h(1-h)) = dlog h - dlog(1-h) :
  S(Li_{3,1}) = S(Li_{2,1}) x [x] + [1-x][x][x] x [1-y] + [1-xy]([x]+[y])^2 x ([y]-[1-y])
  S(Li_{2,1}) = S(Li_{1,1}) x [x] + [1-x][x] x [1-y] + [1-xy]([x]+[y]) x ([y]-[1-y])
  S(Li_{1,1}) = [1-xy][1-x] + [1-x][1-y] + [1-xy]([y]-[1-y])
  S(Li_{2,2}) = S(Li_{1,2}) x [x] + S(Li_{2,1}) x [y]
  S(Li_{1,2}) = [1-xy]([x]+[y]) x [1-x] + S(Li_{1,1}) x [y]
Validated: y=1 -> H_{0,0,1,1}, H_{0,1,0,1}; y=-1 -> -H_{0,0,1,-1}, -H_{0,1,0,-1};
stuffle Li_{1,1}(x,y)+Li_{1,1}(y,x) = Li_1(x)Li_1(y) - Li_2(xy).
"""
from fractions import Fraction as Fr

def ladd(u, v): return tuple(a + b for a, b in zip(u, v))
def lneg(u): return tuple(-a for a in u)

# Laurent polynomials: dict e -> Fr
def padd(p, q, s=1):
    out = dict(p)
    for e, v in q.items():
        out[e] = out.get(e, 0) + s * v
    return {e: v for e, v in out.items() if v}
def pmul(p, q):
    out = {}
    for e1, v1 in p.items():
        for e2, v2 in q.items():
            out[e1 + e2] = out.get(e1 + e2, 0) + v1 * v2
    return {e: v for e, v in out.items() if v}
ONE = {0: Fr(1)}

def sadd(*S):
    out = {}
    for s in S:
        for w, cf in s.items():
            out[w] = padd(out.get(w, {}), cf)
    return {w: cf for w, cf in out.items() if cf}
def sscale(s, q):
    """q: Laurent poly"""
    return {w: pmul(cf, q) for w, cf in s.items()}
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
                out[w] = padd(out.get(w, {}), pmul(cu, cv))
    return {w: cf for w, cf in out.items() if cf}

# ----- basic symbols -----
def li_sym(n, rv, omrv):
    """S(Li_n(g)) = -[1-g][g]^(n-1); rv=[g], omrv=[1-g] as letter vectors."""
    return {(lneg(omrv),) + (rv,) * (n - 1): dict(ONE)}
def log_sym(gv):
    """S(log g) = [g]"""
    return {(gv,): dict(ONE)}

def _t(word, sign=1):
    return {tuple(word): {0: Fr(sign)}}

def li11_sym(X, Y, P, Q, R):
    S = _t((R, P)) ; S = sadd(S, _t((P, Q)), _t((R, Y)), _t((R, Q), -1))
    return S
def li21_sym(X, Y, P, Q, R):
    S11 = li11_sym(X, Y, P, Q, R)
    S = {w + (X,): cf for w, cf in S11.items()}
    S = sadd(S, _t((P, X, Q)), _t((R, X, Y)), _t((R, Y, Y)), _t((R, X, Q), -1), _t((R, Y, Q), -1))
    return S
def li31_sym(X, Y, P, Q, R):
    S21 = li21_sym(X, Y, P, Q, R)
    S = {w + (X,): cf for w, cf in S21.items()}
    S = sadd(S, _t((P, X, X, Q)))
    XY2 = [((R, X, X), 1), ((R, X, Y), 1), ((R, Y, X), 1), ((R, Y, Y), 1)]
    for w, s in XY2:
        S = sadd(S, _t(w + (Y,), s), _t(w + (Q,), -s))
    return S
def li12_sym(X, Y, P, Q, R):
    S11 = li11_sym(X, Y, P, Q, R)
    S = _t((R, X, P)) if False else {}
    S = sadd(S, _t((R, X, P)), _t((R, Y, P)))
    S = sadd(S, {w + (Y,): cf for w, cf in S11.items()})
    return S
def li22_sym(X, Y, P, Q, R):
    S12 = li12_sym(X, Y, P, Q, R)
    S21 = li21_sym(X, Y, P, Q, R)
    S = {w + (X,): cf for w, cf in S12.items()}
    S = sadd(S, {w + (Y,): cf for w, cf in S21.items()})
    return S

# ----- wedges -----
def wedge(u, v, dim):
    out = {}
    for i in range(dim):
        for j in range(i + 1, dim):
            c = u[i] * v[j] - u[j] * v[i]
            if c:
                out[(i, j)] = out.get((i, j), 0) + c
    return out
def delta22(S, dim):
    """True Goncharov delta_{2,2} -> Lambda^2 B_2 (COR4 semantics):
    (a,b,c,d) -> (a^b)^(c^d) with second antisymmetrisation; coefficients Laurent polys.
    Returns dict ((p1,p2)) -> Laurent poly."""
    out = {}
    for w, cf in S.items():
        a, b, cc, d = w
        w1 = wedge(a, cc, dim) if False else wedge(a, b, dim)
        w2 = wedge(cc, d, dim)
        for p1, c1 in w1.items():
            for p2, c2 in w2.items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                out[key] = padd(out.get(key, {}), {e: v * c1 * c2 * sgn for e, v in cf.items()})
    return {k: v for k, v in out.items() if v}

# ----- substitution of letters -----
def sub_letters(S, images, dim):
    """images: dict old-letter-index -> new-letter-vector (length dim).
    Every letter in every word must be a basis vector e_i (or a sum);
    we expand multilinearly."""
    out = {}
    for w, cf in S.items():
        terms = [((), dict(cf))]
        for lett in w:
            new_terms = []
            # decompose lett into basis indices
            opts = [(j, lett[j]) for j in range(len(lett)) if lett[j]]
            for pref, c0 in terms:
                for j, mult in opts:
                    new_terms.append((pref + (images[j],), {e: v * mult for e, v in c0.items()}))
            terms = new_terms
        for pref, c0 in terms:
            out[pref] = padd(out.get(pref, {}), c0)
    return {w: cf for w, cf in out.items() if cf}
