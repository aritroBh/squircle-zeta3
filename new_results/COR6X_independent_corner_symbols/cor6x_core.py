
# cor6x_core.py -- COR6X referee: shared machinery
# Regularized symbol of the NUMBER G(a_1..a_n; x) = I(0; a_1..a_n; x),
# Goncharov recursion derived from scratch:
#   S(I) = sum_i S(I(omit a_i)) (x) ( [prev_i - a_i] - [next_i - a_i] )
#   cyclic convention: prev_1 = endpoint, next_n = startpoint.
# Regularization: a_i -> a_i + eps^i; collisions give letter EPS, killed at the end.

from functools import lru_cache
import sympy as sp
import itertools

EPS = sp.Symbol('EPSREG')   # sentinel, never a true letter

def canon(e):
    """canonical representative of the letter class [e] = [-e] (torsion -1 killed)."""
    e = sp.simplify(e)
    if e is EPS:
        return e
    ne = sp.simplify(-e)
    ke, kn = (str(e).startswith('-'), str(e)), (str(ne).startswith('-'), str(ne))
    return e if ke <= kn else ne

def is_zero_class(e):
    """letter e is zero in V: constant root of unity."""
    e = sp.simplify(e)
    if e is EPS or e.free_symbols:
        return False
    if e in (1, -1, sp.I, -sp.I):
        return True
    try:
        return sp.simplify(e**2 - 1) == 0
    except Exception:
        return False

def dlog(P, Q):
    """log-class of P - Q. P,Q: ('c', expr) endpoint or ('m', expr, i) letter a_i+eps^i."""
    d = sp.simplify(P[1] - Q[1])
    if d == 0:
        if P[0] == 'c' and Q[0] == 'c':
            raise AssertionError("both endpoints equal: integral identically zero")
        return EPS
    return canon(d)

@lru_cache(maxsize=None)
def sym_I(a0, letters, aN):
    """Symbol of regularized I(a0; letters; aN); a0,aN sympy exprs, letters tuple."""
    n = len(letters)
    A0 = ('c', a0); AN = ('c', aN)
    LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
    if n == 0:
        return {(): 1}
    out = {}
    if n == 1:
        e1 = dlog(AN, LL[0]); e0 = dlog(A0, LL[0])
        if not (e1 is EPS or is_zero_class(e1)):
            out[(e1,)] = out.get((e1,), 0) + 1
        if not (e0 is EPS or is_zero_class(e0)):
            out[(e0,)] = out.get((e0,), 0) - 1
        return {k: v for k, v in out.items() if v}
    for i in range(n):
        sub = tuple(l[1] for l in LL[:i] + LL[i+1:])
        Ssub = sym_I(a0, sub, aN)
        prev = AN if i == 0 else LL[i-1]
        nxt = A0 if i == n - 1 else LL[i+1]
        e_prev = dlog(prev, LL[i])
        e_next = dlog(LL[i], nxt)
        for k, v in Ssub.items():
            if not (e_prev is EPS or is_zero_class(e_prev)):
                kk = k + (e_prev,)
                out[kk] = out.get(kk, 0) + v
            if not (e_next is EPS or is_zero_class(e_next)):
                kk = k + (e_next,)
                out[kk] = out.get(kk, 0) - v
    out = {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}
    return out

def sym_G(letters, pt):
    return sym_I(sp.Integer(0), tuple(sp.sympify(a) for a in letters), sp.sympify(pt))

def sym_H(word, pt):
    """Symbol of H_word(pt); Remiddi-Vermaseren H = (-1)^{#ones} G (checked numerically)."""
    sgn = (-1) ** sum(1 for a in word if a == 1)
    return {k: sgn * v for k, v in sym_G(word, pt).items()}

# ---------- symbol algebra ----------
def shuffle_words(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + r for r in shuffle_words(u[1:], v)] + \
           [(v[0],) + r for r in shuffle_words(u, v[1:])]

def sym_shuffle(S, T):
    out = {}
    for m1, c1 in S.items():
        for m2, c2 in T.items():
            for w in shuffle_words(m1, m2):
                out[w] = out.get(w, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}

def sym_add(*Ss):
    out = {}
    for S in Ss:
        for k, v in S.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def sym_scale(S, c):
    return {k: v * c for k, v in S.items() if v * c != 0}

def sym_sub(S, T):
    return sym_add(S, sym_scale(T, -1))

# ---------- wedge / delta ----------
def _wedge(a, b):
    if a == b: return None, 1
    if str(a) < str(b): return (a, b), 1
    return (b, a), -1

def delta11(S):
    """weight-2 symbol -> Lambda^2 V canonical dict."""
    out = {}
    for (a, b), v in S.items():
        w, s = _wedge(a, b)
        if w is None: continue
        out[w] = out.get(w, 0) + s * v
    return {k: v for k, v in out.items() if v}

def delta22(S):
    """weight-4 symbol -> Lambda^2(Lambda^2 V) canonical dict."""
    out = {}
    for (a, b, c, d), v in S.items():
        w1, s1 = _wedge(a, b); w2, s2 = _wedge(c, d)
        if w1 is None or w2 is None or w1 == w2: continue
        s = s1 * s2
        if str(w1) < str(w2):
            key = (w1, w2)
        else:
            key = (w2, w1); s = -s
        out[key] = out.get(key, 0) + s * v
    return {k: v for k, v in out.items() if v}

def _dlog_wedge(a, b, var):
    """coefficient of dlog(a) ^ dlog(b) as a 2-form in the given variables:
       returns (dvar_a ^ dvar_b) = d(a)/a ^ d(b)/b, as a sympy expr per pair of vars."""
    return sp.simplify(sp.diff(a, var[0]) / a * sp.diff(b, var[1]) / b
                       - sp.diff(a, var[1]) / a * sp.diff(b, var[0]) / b)

def integrable(S, vars=None):
    """symbol integrability: for every adjacent pair position, sum c * dlog(wi)^dlog(wj) (x) rest = 0
       as a differential 2-form tensored with the remaining entries."""
    if vars is None:
        fs = set()
        for k in S:
            for t in k: fs |= t.free_symbols
        vars = sorted(fs, key=str)
    bad = {}
    for k, v in S.items():
        for i in range(len(k) - 1):
            for p in range(len(vars)):
                for q in range(p + 1, len(vars)):
                    cf = _dlog_wedge(k[i], k[i+1], (vars[p], vars[q]))
                    if cf != 0:
                        key = (i, p, q, k[:i] + k[i+2:])
                        bad[key] = bad.get(key, 0) + v * cf
    return {k: sp.simplify(v) for k, v in bad.items() if sp.simplify(v) != 0}


# ---------- letter expansion into a prime basis (rational-function letters) ----------
def expand_letter(e):
    """e: sympy rational expression over Q. Returns {irreducible factor (canon. sign): exp}.
       Torsion (+-1) dropped."""
    num, den = sp.fraction(sp.cancel(e))
    out = {}
    for poly_expr, sgn in ((num, 1), (den, -1)):
        if not poly_expr.free_symbols:
            r = sp.Rational(poly_expr)
            for b, m in sp.factorint(r.p).items():
                if b != -1:
                    out[sp.Integer(b)] = out.get(sp.Integer(b), 0) + sgn * m
            for b, m in sp.factorint(r.q).items():
                if b != -1:
                    out[sp.Integer(b)] = out.get(sp.Integer(b), 0) - sgn * m
            continue
        content, fac = sp.factor_list(poly_expr)
        r = sp.Rational(content)
        for b, m in sp.factorint(r.p).items():
            if b != -1:
                out[sp.Integer(b)] = out.get(sp.Integer(b), 0) + sgn * m
        for b, m in sp.factorint(r.q).items():
            if b != -1:
                out[sp.Integer(b)] = out.get(sp.Integer(b), 0) - sgn * m
        for f, m in fac:
            f0 = sp.primitive(f.as_expr())[1]
            if sp.simplify(f0 - 1) == 0:
                continue
            if sp.LC(sp.Poly(f0)) < 0:
                f0 = -f0
            out[f0] = out.get(f0, 0) + sgn * m
    return {k: v for k, v in out.items() if v}

def sym_expanded(S):
    """multilinear expansion of every letter into the prime basis ([prod p^e] = sum e [p])."""
    import itertools as _it
    out = {}
    for mono, c in S.items():
        exps = [expand_letter(t) for t in mono]
        if any(len(e) == 0 for e in exps):
            continue
        for combo in _it.product(*[list(e.items()) for e in exps]):
            newmono = tuple(k for k, v in combo)
            mult = c
            for k, v in combo:
                mult *= v
            out[newmono] = out.get(newmono, 0) + mult
    return {k: v for k, v in out.items() if v}
