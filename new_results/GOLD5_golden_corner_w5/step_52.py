lib = r'''"""gold5_lib.py -- GOLD5 shared library: golden W5 points, NF cobracket/quotient machinery.

Reuses G2b's fixed pipeline (step_01/step_05 modules, loaded by importlib; never edited)
and COR6/G2_M4alpha machinery. Conventions identical to G2b (paper sec:M4):
  obs pairing with X = [2], Y = [phi] (basis letters); Q = [X,[X,Y]], T = [Y,[X,Y]], P = [X,Y].
  delta_word_sym: Goncharov (3,2) cobracket at the symbol level (M4c port, as in G2b step_05).
  Li_{3,2}(x,y) = sum_{n>m>0} x^n y^m/(n^3 m^2)  (paper convention; convergent atoms: |x|<=1, |xy|<=1).
  J((0,0,1/x,0,1/(xy));0->1) = +Li_{3,2}(x,y).
NF elements: tuples of Fractions, (a, b) = a + b*sqrt5.
"""
import sys, os, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
GOLD = BASE + "/GOLD5_golden_corner_w5/"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, GOLD)
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); mod.__name__ = name
    spec.loader.exec_module(mod)
    return mod
s01 = load("g2b_step01", G2B + "step_01.py")
from cor6_sym import sym_J, sadd, sscale, sshuffle, rref
from step_13 import W5Point
from step_12 import br
import g2_fit as gf
from step_10 import load_terms

PHI = (Fr(1, 2), Fr(1, 2))          # (1+sqrt5)/2
TWO = (Fr(2), Fr(0))
BASIS = [TWO, PHI]
SQ5 = (Fr(0), Fr(1))                # sqrt5

def nf_str(a):
    if len(a) == 1: return str(a[0])
    if a[1] == 0: return str(a[0])
    if a[0] == 0: return f"{a[1]}*sqrt5"
    return f"{a[0]}{'+' if a[1] > 0 else ''}{a[1]}*sqrt5"

def golden_point(which):
    if which == 'gol1':   # c0 = 1/phi = (-1+sqrt5)/2
        P = W5Point("c=1/phi (golden1)", [-5, 0, 1], (Fr(-1, 2), Fr(1, 2)), basis=BASIS)
    elif which == 'gol3':  # c0 = sqrt5-2 = phi^-3
        P = W5Point("c=sqrt5-2 (golden3)", [-5, 0, 1], (Fr(-2), Fr(1)), basis=BASIS)
    else:
        raise ValueError(which)
    P.S0 = [2]; P.dimV = 2
    return P

def log(msg):
    print(msg, flush=True)
    with open(GOLD + "progress.log", "a") as f:
        f.write(msg + "\n")

# ---------------- NF versions of G2b's step_05 machinery ----------------
def make_S_J(P):
    F = P.F
    def S_J(word, x, a0=0):
        return sym_J(tuple(word), x, P.vec, a0=a0,
                     sub=F.sub, iszero=lambda a: all(t == 0 for t in a))
    return S_J

def lam2(P, word, x, a0=0):
    """Lambda^2 coordinate ([2]^[phi] coeff) of a 2-letter J-symbol."""
    S = make_S_J(P)(word, x, a0)
    out = Fr(0)
    for (u, v), cf in S.items():
        out += cf * (u[0] * v[1] - u[1] * v[0])
    return out

def hpl_sign(w):
    return Fr(-1) if sum(1 for a in w if a == 1) % 2 else Fr(1)

def signed_words(P):
    """the (word, NF coefficient) pairs of Omega_5 at the point, as J-words with HPL signs."""
    F = P.F
    out = []
    for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
        if len(w) != 5: continue
        sgn = Fr(1) if (e + gf.nz(w)) % 2 == 0 else Fr(-1)
        c0e = F.pow(P.c0, e)
        coef1 = F.mul(c0e, F.scale(q * hpl_sign(w), F.one()))
        out.append((w, coef1))
        coef2 = F.mul(F.mul(c0e, F.scale(q * sgn, F.one())), F.scale(hpl_sign(gf.neg(w)), F.one()))
        out.append((gf.neg(w), coef2))
    return out

def delta_word_sym(P, w, e, c):
    """symbol-level (3,2) cobracket of J((e,w;c)) : dict left-3-word (exponent vectors) -> NF element."""
    F = P.F
    S_J = make_S_J(P)
    pts = (e,) + tuple(w) + (c,)
    out = {}
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = lam2(P, gap, ci, ei)
                if not lam: continue
                for lw, cf in S_J(ret, c, e).items():
                    out[lw] = F.add(out.get(lw, F.zero()), F.scale(cf * lam, F.one()))
            else:
                lam = lam2(P, ret, c, e)
                if not lam: continue
                for lw, cf in S_J(gap, ci, ei).items():
                    out[lw] = F.sub(out.get(lw, F.zero()), F.scale(cf * lam, F.one()))
    return {k: v for k, v in out.items() if any(v)}

def delta_minus_5sym(S5):
    """delta^-(S) on the 5-symbol (exponent-vector words, NF coefficients)."""
    out = {}
    for w, cf in S5.items():
        a, b = w[3], w[4]
        lam = a[0] * b[1] - a[1] * b[0]
        if lam:
            lw = w[:3]
            out[lw] = tuple(x + lam * y for x, y in zip(out.get(lw, (Fr(0),) * len(cf)), cf))
        a, b = w[0], w[1]
        lam = a[0] * b[1] - a[1] * b[0]
        if lam:
            lw = w[2:]
            out[lw] = tuple(x - lam * y for x, y in zip(out.get(lw, (Fr(0),) * len(cf)), cf))
    return {k: v for k, v in out.items() if any(v)}

# ---- quotient T^3/D^3 over 2 letters with NF coefficients ----
def t3_basis():
    from itertools import product
    return list(product(range(2), repeat=3))

def d3_rows():
    from itertools import product
    rows = []
    for a in range(2):
        for bc in product(range(2), repeat=2):
            rows.append(sshuffle({(a,): Fr(1)}, {bc: Fr(1)}))
    return rows

_D3 = None
def d3_rref():
    global _D3
    if _D3 is None:
        B = t3_basis()
        Dmat, piv = rref([[s.get(w, Fr(0)) for w in B] for s in d3_rows()])
        _D3 = (Dmat, piv)
    return _D3

def quotient_coords_nf(SE, dimF=2):
    """SE: {index-word -> NF element} at weight 3 -> coords in T^3/D^3 (free words, NF)."""
    Dmat, piv = d3_rref()
    B = t3_basis()
    free = [j for j in range(len(B)) if j not in piv]
    v = [tuple(SE.get(w, (Fr(0),) * dimF)) for w in B]
    for i, pj in enumerate(piv):
        if any(v[pj]):
            q = v[pj]
            v = [tuple(vv[t] - q[t] * Dmat[i][j] for t in range(dimF)) for j, vv in enumerate(v)]
    for j in range(len(B)):
        if j not in free:
            assert not any(v[j]), "not integrable"
    return tuple(v[j] for j in free), free

def ft_expand_sym_nf(S, dim):
    """S: {exp-vector-word -> NF}. Expand letters to basis indices; NF coeffs multiply Fractions."""
    out = {}
    for w, cf in S.items():
        terms = [((), Fr(1))]
        for l in w:
            terms = [(p + (j,), c * l[j]) for p, c in terms for j in range(dim) if l[j]]
        for p, c in terms:
            nc = tuple(c * x for x in cf)
            cur = out.get(p)
            out[p] = nc if cur is None else tuple(a + b for a, b in zip(cur, nc))
    return {w: c for w, c in out.items() if any(c)}

def li32_word(P, x, y):
    """J-word of Li_{3,2}(x,y): (0,0,1/x,0,1/(xy)), NF elements."""
    F = P.F
    return (F.zero(), F.zero(), F.inv(x), F.zero(), F.inv(F.mul(x, y)))

def li3_sym_nf(P, z):
    """S(Li3(z)) = -[1-z][z][z], exponent-vector letters, NF coeff."""
    F = P.F
    vz = P.vec(z); v1z = P.vec(F.sub(F.one(), z))
    return {(v1z, vz, vz, vz): (Fr(-1),) + (Fr(0),) * (F.d - 1)}

def li32_sym_nf(P, x, y):
    return make_S_J(P)(li32_word(P, x, y), P.F.one(), P.F.zero())

def obs_of_SE(SE):
    """pairing of an expanded weight-5 symbol {index-word -> Fr or NF}: returns [obs_Q, obs_T] (NF)."""
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    P = br(X, Y); Q = br(X, br(X, Y)); T = br(Y, br(X, Y))
    out = []
    for A in (Q, T):
        tot = None
        for u, cu in A.items():
            for v, dv in P.items():
                for word, sgn in ((u + v, Fr(1)), (v + u, Fr(-1))):
                    c = SE.get(word)
                    if c is None: continue
                    cf = c if isinstance(c, tuple) else (c,)
                    cc = tuple(sgn * cu * dv * x for x in cf)
                    tot = cc if tot is None else tuple(a + b for a, b in zip(tot, cc))
        out.append(tot if tot is not None else ())
    return out
'''
open('gold5_lib.py', 'w').write(lib)
import sys
sys.path.insert(0, '.')
import gold5_lib as G
from fractions import Fraction as Fr
print("lib loads OK")
qc, free = G.quotient_coords_nf({(0,1,1): (Fr(3), Fr(-2))})
print("quotient test:", [[str(x) for x in q] for q in qc], "free:", free)
