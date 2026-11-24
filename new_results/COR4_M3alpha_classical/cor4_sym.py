"""cor4_sym.py -- exact symbol algebra for weight-4 HPLs over letters {0,1,-1}
at arguments +/-c, with rational-function coefficients c^e tracked as scalars.

Conventions (checked against S(Li_n(x)) = -[1-x] [x]^(n-1)):
  iterated integral I(0; b_1,...,b_n; x) = int_{0<t_1<...<t_n<x} w_{b_1}(t_1)...w_{b_n}(t_n),
  w_0 = dlog t, w_1 = dlog(1/(1-t)), w_{-1} = dlog(1+t).
  H_{a_1,...,a_n}(x) = I(0; a_n,...,a_1; x), so
  S(H_{a_1,...,a_n}(x)) = g(a_n) x ... x g(a_1), g(0)=[x], g(1)=-[1-x], g(-1)=[1+x].

A letter is a tuple of Fractions = exponent vector over a fixed multiplicative basis.
A symbol is dict: word(tuple of letters) -> Fraction coefficient.
A "function symbol" tracks the scalar c^e: dict (e, word) -> Fraction.
"""
from fractions import Fraction as Fr

# ---------- letter vectors ----------
def ladd(u, v): return tuple(a + b for a, b in zip(u, v))
def lsub(u, v): return tuple(a - b for a, b in zip(u, v))
def lneg(u): return tuple(-a for a in u)
def lscale(s, u): return tuple(s * a for a in u)

# ---------- symbols (dict word -> Fr) ----------
def sadd(*S):
    out = {}
    for s in S:
        for w, c in s.items():
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def sscale(s, c):
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

# ---------- HPL / Li symbols ----------
def hpl_sym(word, g):
    """S(H_word(x)) with g = {0:[x], 1:-[1-x], -1:[1+x]} mapped into letter vectors."""
    return {tuple(g[a] for a in reversed(word)): Fr(1)}

def li_sym(n, rv, omrv):
    """S(Li_n(r)) = -[1-r] [r]^(n-1).  rv=[r], omrv=[1-r] as letter vectors."""
    return {(lneg(omrv),) + (rv,) * (n - 1): Fr(1)}

# ---------- wedges and the (2,2) cobracket ----------
def wedge(u, v, dim):
    """u wedge v as dict (i,j), i<j -> coefficient."""
    out = {}
    for i in range(dim):
        for j in range(i + 1, dim):
            c = u[i] * v[j] - u[j] * v[i]
            if c:
                out[(i, j)] = out.get((i, j), 0) + c
    return out

def delta22(S, dim):
    """Antisymmetrised (2,2) component: (a,b,c,d) -> (a^b) x (c^d).
    Kills shuffle products and classical Li4.  Returns dict ((p1,p2)) -> Fr."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        w1 = wedge(a, b, dim)
        w2 = wedge(c, d, dim)
        for p1, c1 in w1.items():
            for p2, c2 in w2.items():
                out[(p1, p2)] = out.get((p1, p2), 0) + cf * c1 * c2
    return {k: v for k, v in out.items() if v}

def delta22_func(FS, dim):
    """Same for a function symbol dict (e, word) -> Fr.
    Returns dict wedgepair -> dict e -> Fr (Laurent polynomial coefficients)."""
    out = {}
    for (e, w), cf in FS.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                d = out.setdefault((p1, p2), {})
                d[e] = d.get(e, 0) + cf * c1 * c2
    return {k: {e: v for e, v in d.items() if v} for k, d in out.items()
            if any(d.values())}

# ---------- function symbols: dict (e, word) -> Fr ----------
def fadd(*S):
    out = {}
    for s in S:
        for k, c in s.items():
            out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}

def fscale_c(FS, k):
    """Multiply a function symbol by the scalar c^k."""
    return {(e + k, w): c for (e, w), c in FS.items()}

def fscale(FS, q):
    return {k: v * q for k, v in FS.items() if v * q}


# ---------- specialisation ----------
def make_letter_map(images):
    """images: list of new-basis vectors = images of the old basis letters.
    Returns function mapping an old letter vector to a new letter vector."""
    def f(u):
        dim = len(images[0])
        return tuple(sum(u[j] * images[j][k] for j in range(len(u))) for k in range(dim))
    return f

def specialize_func(FS, lmap, cpow):
    """FS: function symbol dict (e, word) -> Fr.
    lmap: letter vector -> new letter vector.
    cpow: dict e -> exact field element for c0^e (Fr or (Fr,Fr) pair).
    Returns dict word -> field-element coefficient (same type as cpow values)."""
    out = {}
    for (e, w), cf in FS.items():
        nw = tuple(lmap(v) for v in w)
        c0e = cpow[e]
        if isinstance(c0e, tuple):
            val = (cf * c0e[0], cf * c0e[1])
        else:
            val = cf * c0e
        out[nw] = out.get(nw, 0) + val if not isinstance(val, tuple) else (
            (out.get(nw, (0, 0))[0] + val[0], out.get(nw, (0, 0))[1] + val[1]))
    if out and isinstance(next(iter(out.values())), tuple):
        return {w: c for w, c in out.items() if c != (Fr(0), Fr(0))}
    return {w: c for w, c in out.items() if c != 0}

def delta22_spec(S, dim):
    """(2,2) cobracket for specialised symbols whose coefficients may be (Fr,Fr) pairs."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                key = (p1, p2)
                if isinstance(cf, tuple):
                    cur = out.get(key, (Fr(0), Fr(0)))
                    out[key] = (cur[0] + cf[0] * c1 * c2, cur[1] + cf[1] * c1 * c2)
                else:
                    out[key] = out.get(key, 0) + cf * c1 * c2
    if out and isinstance(next(iter(out.values())), tuple):
        return {k: v for k, v in out.items() if v != (Fr(0), Fr(0))}
    return {k: v for k, v in out.items() if v}

def canon(S, dim):
    """Expand a symbol whose letters are arbitrary vectors into the tensor
    power of the basis (e_0,...,e_{dim-1}); faithful normal form of V^{x n}."""
    out = {}
    for w, cf in S.items():
        terms = [((), cf)]
        for v in w:
            terms = [(pref + (j,), c * v[j])
                     for pref, c in terms for j in range(dim) if v[j]]
        for pref, c in terms:
            out[pref] = out.get(pref, 0) + c
    return {k: v for k, v in out.items() if v}


# ---------- corrected (2,2) cobracket: target Lambda^2(B_2) ----------
def delta22L2(S, dim):
    """Goncharov (2,2) component into Lambda^2(B_2):
    (l1,l2,l3,l4) -> (l1^l2)^(l3^l4)  (antisymmetrised within AND between factors).
    Returns dict ((ij),(kl)) with (ij)<(kl) -> Fr.  Kills Li4 and all shuffle products."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                out[key] = out.get(key, 0) + cf * c1 * c2 * sgn
    return {k: v for k, v in out.items() if v}

def delta22L2_func(FS, dim):
    """Function-level version: dict wedgepair -> dict e -> Fr."""
    out = {}
    for (e, w), cf in FS.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                dd = out.setdefault(key, {})
                dd[e] = dd.get(e, 0) + cf * c1 * c2 * sgn
    return {k: {e: v for e, v in dd.items() if v} for k, dd in out.items() if any(dd.values())}

def delta22L2_spec(S, dim):
    """Specialised version; coefficients may be Fr or (Fr,Fr)."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                v = c1 * c2 * sgn
                if isinstance(cf, tuple):
                    cur = out.get(key, (Fr(0), Fr(0)))
                    out[key] = (cur[0] + cf[0]*v, cur[1] + cf[1]*v)
                else:
                    out[key] = out.get(key, 0) + cf * v
    if out and isinstance(next(iter(out.values())), tuple):
        return {k: v for k, v in out.items() if v != (Fr(0), Fr(0))}
    return {k: v for k, v in out.items() if v}

# ---------- (3,1) and (1,3) components via Lie-word pairings ----------
def lie_bracket(words):
    """[[l1,l2],l3] as dict word->coeff from a tuple of 3 letters."""
    l1, l2, l3 = words
    out = {}
    for w, c in [((l1,l2,l3),1),((l2,l1,l3),-1),((l3,l1,l2),-1),((l3,l2,l1),1)]:
        out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def delta31(S, dim):
    """Pair first three letters with the Lie bracket [[x,y],z], tensor the fourth.
    Returns dict (basis3word, letterindex) -> Fr where basis3word is a word in the
    letter-basis indices after expanding."""
    out = {}
    for w, cf in S.items():
        l1, l2, l3, l4 = w
        for bw, bc in lie_bracket((l1, l2, l3)).items():
            # expand all four letters into basis indices
            for j4 in range(dim):
                if not l4[j4]:
                    continue
                terms = [((), cf * bc * l4[j4])]
                for v in bw:
                    terms = [(p + (j,), c * v[j]) for p, c in terms for j in range(dim) if v[j]]
                for p, c in terms:
                    out[(p, j4)] = out.get((p, j4), 0) + c
    return {k: v for k, v in out.items() if v}

def delta13(S, dim):
    out = {}
    for w, cf in S.items():
        l1, l2, l3, l4 = w
        for bw, bc in lie_bracket((l2, l3, l4)).items():
            for j1 in range(dim):
                if not l1[j1]:
                    continue
                terms = [((), cf * bc * l1[j1])]
                for v in bw:
                    terms = [(p + (j,), c * v[j]) for p, c in terms for j in range(dim) if v[j]]
                for p, c in terms:
                    out[(j1, p)] = out.get((j1, p), 0) + c
    return {k: v for k, v in out.items() if v}

# ---------- Q(sqrt2) S-unit utilities ----------
def _cmul(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def _cdiv(a, b):
    p, q = b; n = p*p - 2*q*q
    return _cmul(a, (p/n, -q/n))
_SQ2 = (Fr(0), Fr(1)); _A1 = (Fr(1), Fr(1))

def sunit_vec_sq2(z):
    """z=(p,q)=p+q sqrt2 in Q(sqrt2)^*.  If |z| = 2^{j/2} (1+sqrt2)^m with j,m in Z,
    return (Fr(j,2), Fr(m)) [sign/torsion ignored]; else None. Exact verification."""
    if z == (Fr(0), Fr(0)): return None
    p, q = z
    N = p*p - 2*q*q
    if N == 0: return None
    s = 1 if N > 0 else -1
    j2 = None
    n = N if N > 0 else -N
    # n must be a power of 2 (as a rational): numerator & denominator powers of 2
    num, den = n.numerator, n.denominator
    if num & (num - 1) or den & (den - 1):
        return None
    import math
    j = num.bit_length() - den.bit_length()  # n = 2^j
    # z = 2^{j/2} u : compute u = z / 2^{j/2}
    if j % 2 == 0:
        f2 = (Fr(2) ** (j // 2), Fr(0))
    else:
        f2 = _cmul(_SQ2, (Fr(2) ** ((j - 1) // 2), Fr(0)))
    u = _cdiv(z, f2)
    # u must be +- (1+sqrt2)^m
    import math
    uf = float(u[0]) + float(u[1]) * math.sqrt(2)
    if uf == 0: return None
    m = int(round(math.log(abs(uf)) / math.log(1 + math.sqrt(2))))
    uu = (Fr(1), Fr(0))
    for _ in range(abs(m)):
        uu = _cmul(uu, _A1) if m > 0 else _cdiv(uu, _A1)
    if u != uu and u != (-uu[0], -uu[1]):
        return None
    return (Fr(j, 2), Fr(m))
