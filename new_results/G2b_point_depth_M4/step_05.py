"""step_05.py -- G2b Task 4 (exact part): the (3,2) cobracket of R4^{(5)}(1/2) at the SYMBOL level
and the Li_{3,2} correction fixed by it.

Conventions (X = [2], Y = [3], basis-index words over {0,1}):
  delta_word (M4c w5_cobracket.py, ported): for a 5-letter J-word with points (e,w,c),
    Delta = sum_{k=2,i} S(J(retained; e->c)) * Lam(gap_i)  -  sum_{k=3,i} S(J(gap_i; ei->ci)) * Lam(ret_i)
  with Lam = the Lambda^2 coordinate [2]^ [3] of the 2-letter sub-integral's symbol (exact).
  The left 3-letter symbols are exact exponent-vector tensors (COR6 sym_J, a0=ei, x=ci).

EXACT CHECK 1: the assembled Delta_sym equals delta^-(Omega_5(1/2)) computed directly from the
5-symbol, delta^-(S) = sum_w s(w)[ w[:3]*Lam(w[3:]) - w[2:]*Lam(w[:2]) ].  (Validates the port.)

Then Delta_sym = sum_z c_z * S(Li3(z))  mod shuffle products  (2-dim quotient T^3/D^3):
the c_z are solved EXACTLY.  Same for every Li32 atom.  N5 = sum d_j Li32(x_j,y_j) is solved
so that sum d_j c(atom_j) = c(D) exactly (rational); obs(N5) == (41,-33/2) is re-verified.
"""
import sys, os, json, time
from fractions import Fraction as Fr
from math import gcd
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
_spec = importlib.util.spec_from_file_location("g2b_step01", G2B + "step_01.py")
_s01 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_s01)
ft_expand_sym = _s01.ft_expand_sym; log = _s01.log
from cor6_sym import sym_J, sadd, sscale, sshuffle
from step_13 import W5Point
import g2_fit as gf
from step_10 import load_terms

t0 = time.time()
C0 = Fr(1, 2)

# ---- exponent vectors over {2,3} for rationals (exact)
def vec_q(el):
    if el == 0: raise ValueError("zero letter")
    el = abs(el)
    e2, e3 = 0, 0
    n, d = el.numerator, el.denominator
    while n % 2 == 0: n //= 2; e2 += 1
    while d % 2 == 0: d //= 2; e2 -= 1
    while n % 3 == 0: n //= 3; e3 += 1
    while d % 3 == 0: d //= 3; e3 -= 1
    if n != 1 or d != 1: raise ValueError(f"not a {{2,3}}-unit: {el}")
    return (Fr(e2), Fr(e3))

def S_J(word, x, a0=0):
    return sym_J(tuple(word), x, vec_q, a0=a0,
                 sub=lambda a, b: a - b, iszero=lambda a: a == 0)

def lam2(word, x, a0=0):
    """Lambda^2 coordinate (coefficient of [2]^[3]) of a 2-letter J-integral's symbol."""
    S = S_J(word, x, a0)
    out = Fr(0)
    for (u, v), cf in S.items():
        out += cf * (u[0] * v[1] - u[1] * v[0])
    return out

def hpl_sign(w):
    return Fr(-1) if sum(1 for a in w if a == 1) % 2 else Fr(1)

def signed_words():
    """the 146 (word, coefficient) pairs of Omega_5(1/2), as J-words with HPL signs."""
    out = []
    for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
        if len(w) != 5: continue
        sgn = Fr(1) if (e + gf.nz(w)) % 2 == 0 else Fr(-1)
        c0e = C0 ** e
        out.append((w, q * c0e * hpl_sign(w)))
        out.append((gf.neg(w), q * c0e * sgn * hpl_sign(gf.neg(w))))
    return out

def delta_word_sym(w, e, c):
    """symbol-level cobracket: dict left-3-word (exponent vectors) -> rational (Lambda^2 coord)."""
    pts = (e,) + tuple(w) + (c,)
    out = {}
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = lam2(gap, ci, ei)
                if not lam: continue
                for lw, cf in S_J(ret, c, e).items():
                    out[lw] = out.get(lw, Fr(0)) + cf * lam
            else:
                lam = lam2(ret, c, e)
                if not lam: continue
                for lw, cf in S_J(gap, ci, ei).items():
                    out[lw] = out.get(lw, Fr(0)) - cf * lam
    return {k: v for k, v in out.items() if v}

def omega5_half():
    P = W5Point("c=1/2", [-1, 1], (Fr(1, 2),))
    return P.omega5()

def delta_minus_5sym(S5):
    """delta^-(S) = sum_w s(w)[w[:3]*Lam(w[3:]) - w[:2]*Lam(w[:2])] on the 5-symbol (exponent vecs)."""
    out = {}
    for w, cf in S5.items():
        a, b = w[3], w[4]
        lam = a[0] * b[1] - a[1] * b[0]
        cfq = cf[0] if isinstance(cf, tuple) else cf   # F = Q: NF coefficient is a 1-tuple
        if lam:
            lw = w[:3]
            out[lw] = out.get(lw, Fr(0)) + cfq * lam
        a, b = w[0], w[1]
        lam = a[0] * b[1] - a[1] * b[0]
        if lam:
            lw = w[2:]
            out[lw] = out.get(lw, Fr(0)) - cfq * lam
    return {k: v for k, v in out.items() if v}

# ---- quotient T^3/D^3 over 2 letters (exact)
def t3_basis():
    from itertools import product
    return list(product(range(2), repeat=3))

def d3_rows():
    """shuffle generators at weight 3 over 2 letters, expanded."""
    from itertools import product
    rows = []
    for a in range(2):
        for bc in product(range(2), repeat=2):
            sh = sshuffle({(a,): Fr(1)}, {bc: Fr(1)})
            rows.append(sh)
    return rows

def quotient_coords(SE):
    """SE: {index-word -> Fr} at weight 3 -> coords in T^3/D^3 (2-dim), via rref of D^3."""
    SE = {w: (cf[0] if isinstance(cf, tuple) else cf) for w, cf in SE.items()}
    from cor6_sym import rref
    B = t3_basis()
    def vecof(S): return [S.get(w, Fr(0)) for w in B]
    Dmat, piv = rref([vecof(s) for s in d3_rows()])
    free = [j for j in range(len(B)) if j not in piv]
    v = vecof(SE)
    for i, pj in enumerate(piv):
        if v[pj]:
            q = v[pj]
            v = [a - q * b for a, b in zip(v, Dmat[i])]
    assert all(v[j] == 0 for j in range(len(B)) if j not in free), "not integrable (not in D^3 span?)"
    return tuple(v[j] for j in free), free

def li3_sym(z):
    """S(Li3(z)) = -[1-z][z][z], expanded to basis-index words."""
    S = {(vec_q(1 - z), vec_q(z), vec_q(z)): (Fr(-1),)}
    return ft_expand_sym(S, 2)

def li32_word(x, y):
    return (0, 0, Fr(1) / x, 0, Fr(1) / (x * y))

def li32_sym(x, y):
    """symbol of Li_{3,2}(x,y) = J((0,0,1/x,0,1/(xy)); 0->1), up to the global sign (see note)."""
    return S_J(li32_word(x, y), 1, 0)

def cobracket_obs_from_cvec(cvec_by_z, zbasis):
    return None

if __name__ == '__main__':
    log(f"step_05 start {time.strftime('%H:%M:%S')}")
    # ---- assemble Delta_sym from the 146 signed words
    Delta = {}
    for w, coef in signed_words():
        dw = delta_word_sym(w, 0, C0)
        for lw, v in dw.items():
            Delta[lw] = Delta.get(lw, Fr(0)) + coef * v
    Delta = {k: v for k, v in Delta.items() if v}
    log(f"Delta_sym assembled: {len(Delta)} left words ({time.time()-t0:.0f}s)")
    # ---- EXACT CHECK 1: delta_word == -delta^-(Omega_5(1/2)) in the quotient T^3/D^3.
    # (The symbol is a coalgebra map after reversing tensor order: verified per word at debug time;
    # here we check the total: Delta + delta^-(5-symbol) must be a shuffle product.)
    S5 = omega5_half()
    D5 = delta_minus_5sym(S5)
    keys = set(Delta) | set(D5)
    diff = {k: Delta.get(k, Fr(0)) + D5.get(k, Fr(0)) for k in keys}
    diff = {k: v for k, v in diff.items() if v}
    e0, e1 = (Fr(1), Fr(0)), (Fr(0), Fr(1))
    from itertools import product as _ip
    B3 = list(_ip([e0, e1], repeat=3))
    from cor6_sym import rref
    Drows = []
    for a in (e0, e1):
        for bc in _ip([e0, e1], repeat=2):
            sh = sshuffle({(a,): Fr(1)}, {bc: Fr(1)})
            Drows.append([sh.get(w, Fr(0)) for w in B3])
    Dmat, piv = rref([r[:] for r in Drows])
    v = [diff.get(w, Fr(0)) for w in B3]
    M2, _ = rref([r[:] for r in Drows] + [v])
    chk1 = len(M2) == len(Dmat)
    log(f"EXACT CHECK 1 (delta_word == -delta^-(5-sym) mod shuffles, quotient T^3/D^3): "
        f"{'PASS' if chk1 else 'FAIL'}")
    # ---- quotient coords
    DeltaE = ft_expand_sym({w: (cf,) for w, cf in Delta.items()}, 2)
    qc, free = quotient_coords(DeltaE)
    log(f"Delta Lie_3 quotient coords (free words {[t3_basis()[j] for j in free]}): {[str(q) for q in qc]}")
    # ---- Li3 basis atoms: solve Delta = c1 S(Li3(1/3)) + c2 S(Li3(-1/2)) in the quotient
    zcands = [Fr(1, 3), Fr(-1, 2), Fr(2, 3), Fr(-1, 3), Fr(1, 4), Fr(3, 4), Fr(1, 9), Fr(-1, 8), Fr(-2), Fr(-1)]
    zvecs = {}
    for z in zcands:
        zvecs[z] = quotient_coords(li3_sym(z))[0]
    log("Li3 quotient coords: " + str({str(z): [str(q) for q in v] for z, v in zvecs.items()}))
    # pick two independent
    import itertools
    zpair = None
    for z1, z2 in itertools.combinations(zcands, 2):
        a, b = zvecs[z1], zvecs[z2]
        if a[0] * b[1] - a[1] * b[0]:
            zpair = (z1, z2); break
    z1, z2 = zpair
    a, b = zvecs[z1], zvecs[z2]
    det = a[0] * b[1] - a[1] * b[0]
    c1 = (qc[0] * b[1] - qc[1] * b[0]) / det
    c2 = (a[0] * qc[1] - a[1] * qc[0]) / det
    log(f"Delta's Li3 part: c_{z1} = {c1}, c_{z2} = {c2}   (EXACT)")
    # ---- Li32 atoms: quotient c-vectors and obs pairings
    # atoms: x, y, xy must all be {2,3}-unit-equation solutions (so every letter is smooth)
    sols = [Fr(-8), Fr(-3), Fr(-2), Fr(-1), Fr(-1, 2), Fr(-1, 3), Fr(-1, 8), Fr(1, 9), Fr(1, 4),
            Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(8, 9), Fr(9, 8), Fr(4, 3), Fr(3, 2),
            Fr(2), Fr(3), Fr(4), Fr(9), Fr(1)]
    atoms = []
    for x in sols:
        for y in sols:
            if x == 1 or y == 1 and y != 1: continue
            if x * y in sols and abs(x * y) <= 1 and x != 1:
                atoms.append((x, y))
    atoms = sorted(set(atoms))
    from step_12 import br
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    P = br(X, Y); Q = br(X, br(X, Y)); T = br(Y, br(X, Y))
    def obs_of(SE):
        SE = {w: (cf[0] if isinstance(cf, tuple) else cf) for w, cf in SE.items()}
        out = []
        for A in (Q, T):
            tot = Fr(0)
            for u, cu in A.items():
                for v, dv in P.items():
                    tot += cu * dv * (SE.get(u + v, Fr(0)) - SE.get(v + u, Fr(0)))
            out.append(tot)
        return out
    adata = {}
    for x, y in atoms:
        dw = delta_word_sym(li32_word(x, y), 0, 1)
        dE = ft_expand_sym({w: (cf,) for w, cf in dw.items()}, 2)
        qcw = quotient_coords(dE)[0]
        S_at = ft_expand_sym({w: (cf,) for w, cf in li32_sym(x, y).items()}, 2)
        obs = obs_of(S_at)
        adata[(str(x), str(y))] = {'cvec': [str(q) for q in qcw], 'obs': [str(o) for o in obs]}
        log(f"Li32({x},{y}): cvec = {[str(q) for q in qcw]}, obs = {[str(o) for o in obs]}")
    # ---- solve N5 = sum d_j atom_j with Delta-cobracket == Delta's, in quotient coords
    target = qc
    log(f"solve target (quotient coords): {[str(q) for q in target]} = {c1}*Li3(1/3) + {c2}*Li3(-1/2)")
    best = None
    keys = list(adata.keys())
    for i, j in itertools.combinations(range(len(keys)), 2):
        k1, k2 = keys[i], keys[j]
        a1, a2 = [Fr(s) for s in adata[k1]['cvec']], [Fr(s) for s in adata[k2]['cvec']]
        det = a1[0] * a2[1] - a1[1] * a2[0]
        if not det: continue
        d1 = (target[0] * a2[1] - target[1] * a2[0]) / det
        d2 = (a1[0] * target[1] - a1[1] * target[0]) / det
        h = max(abs(d1.numerator), abs(d1.denominator), abs(d2.numerator), abs(d2.denominator))
        if best is None or h < best[0]:
            best = (h, (k1, d1), (k2, d2))
    log(f"N5 candidate (min-height exact c-solve): {best}")
    # also verify the symbol-level pairing of the candidate
    if best:
        (k1, d1), (k2, d2) = best[1], best[2]
        o1 = [Fr(s) for s in adata[k1]['obs']]; o2 = [Fr(s) for s in adata[k2]['obs']]
        on5 = [d1 * o1[0] + d2 * o2[0], d1 * o1[1] + d2 * o2[1]]
        log(f"obs(N5) = {[str(o) for o in on5]} (target obs(Omega_5) = (41, -33/2)): "
            f"{'MATCH' if on5 == [Fr(41), Fr(-33, 2)] else 'MISMATCH'}")
    json.dump({'Delta_qc': [str(q) for q in qc], 'free': [str(t3_basis()[j]) for j in free],
               'c_li3': {str(z1): str(c1), str(z2): str(c2)},
               'atoms': {f'Li32({k[0]},{k[1]})': v for k, v in adata.items()},
               'N5_candidate': [[best[1][0], str(best[1][1])], [best[2][0], str(best[2][1])]] if best else None},
              open(G2B + 'cobracket_hex.json', 'w'), indent=1)
    log(f"step_05 done ({time.time()-t0:.0f}s)")
