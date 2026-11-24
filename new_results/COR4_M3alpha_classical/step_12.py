cor4_sym = '''"""cor4_sym.py -- exact symbol algebra for weight-4 HPLs over letters {0,1,-1}
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
def specialize_func(FS, letter_map, cpow):
    """letter_map: basis index -> new letter vector (or None to drop dimension).
    cpow: dict e -> exact field element for c0^e (Fr for rational c0).
    Returns dict word -> field-element coefficient."""
    out = {}
    for (e, w), cf in FS.items():
        nw = tuple(letter_map(i, v) for i, v in enumerate(w)) if callable(letter_map) else None
        # letter_map applied componentwise below instead
        raise RuntimeError('use specialize below')
    return out

def spec_word(w, basis_images):
    """basis_images: list where entry j is a function lettervector->newvector."""
    return tuple(basis_images[j](v) for j, v in enumerate(w))
'''
cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
with open(cor + "cor4_sym.py", "w") as f:
    f.write(cor4_sym)
print("written", len(cor4_sym))
