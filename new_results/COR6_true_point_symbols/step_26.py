lib = r'''"""cor6_sym.py -- COR6: TRUE number-level symbols of H_w(c0) including constant letters.

Conventions (control-tested):
  J(a1..an; x): iterated integral, leftmost letter = OUTERMOST:
     J(a1..an;x) = int_0^x dt/(t-a1) J(a2..an; t),  J(empty;x)=1,
     J(0;x) = log x (regularized), J(a;x) = log(1-x/a).
  H_w(x) = (-1)^{# of 1s in w} J(w; x)   (HPL, leftmost = outermost).
  Symbol recursion (Goncharov, shuffle-regularized):
     S(J(a1..an;x)) = sum_j  [l_j] (x) S(J(a1..hat a_j..an;x)),
     l_j = (a_{j+1}-a_j)/(a_{j-1}-a_j),  a_0 = 0, a_{n+1} = x;
     zero numerators/denominators are DROPPED (regularization).
  Tensors: dict word(tuple of letter vectors) -> Fraction.
  Letters are vectors in V = Q^r (exact S-unit exponent vectors).
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

# ---------------- Goncharov symbol recursion ----------------
# field elements are arbitrary hashables; vec(x) returns a letter vector or None for 0.
_MEMO = {}
def sym_J(word, x, vec, a0=0, prepend=True, inv=False):
    """Symbol of J(word; x) with endpoints a0=0 and x.
    vec: element -> letter vector (tuple of Fr) or None if element == 0.
    prepend: letter on the left (q1s convention) else on the right.
    inv: use reciprocal letter ratio."""
    n = len(word)
    key = (word, x, a0, prepend, inv)
    if key in _MEMO: return _MEMO[key]
    if n == 0:
        _MEMO[key] = {}; return {}
    res = {}
    for j in range(n):
        p = a0 if j == 0 else word[j - 1]
        q = x if j == n - 1 else word[j + 1]
        wj = word[j]
        # letter = [(q - wj)/(p - wj)] (as vector); None entries dropped
        num = None if q == wj else (q - wj) if isinstance(q, (int, Fr)) and isinstance(wj, (int, Fr)) else ("sub", q, wj)
        den = None if p == wj else (p - wj) if isinstance(p, (int, Fr)) and isinstance(wj, (int, Fr)) else ("sub", p, wj)
        yield_letter = True
        # delegate vector computation to callback vecdiff
        lv = vecdiff(q, wj, p, vec, inv)
        if lv is None:
            continue  # zero letter -> whole term vanishes
        sub = word[:j] + word[j + 1:]
        Ssub = sym_J(sub, x, vec, a0, prepend, inv)
        if prepend:
            for w2, c2 in Ssub.items():
                w3 = (lv,) + w2
                res[w3] = res.get(w3, 0) + c2
        else:
            for w2, c2 in Ssub.items():
                w3 = w2 + (lv,)
                res[w3] = res.get(w3, 0) + c2
    res = {w: c for w, c in res.items() if c}
    _MEMO[key] = res
    return res

def vecdiff(q, wj, p, vec, inv):
    """letter vector of [(q-wj)/(p-wj)]; None if the letter is 0 (0/0 or empty).
    vec maps a difference token to a vector; uses rational arithmetic when possible."""
    num_zero = (q == wj)
    den_zero = (p == wj)
    if num_zero and den_zero:
        return None
    vn = None if num_zero else vec(("d", q, wj))
    vd = None if den_zero else vec(("d", p, wj))
    if vn is None and vd is None:
        return None
    if inv:
        vn, vd = vd, vn
    if vn is None: return lneg(vd)
    if vd is None: return vn
    return lsub(vn, vd)
'''
with open("cor6_sym.py", "w") as f:
    f.write(lib)
print("written", len(lib))
