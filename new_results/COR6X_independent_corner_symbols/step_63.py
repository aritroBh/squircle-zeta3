t = r'''
# step_08.py -- for each cyclotomic point: letter-group presentation (exact) + delta22 verdict
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp, json, re
import mpmath as mp
mp.mp.dps = 80

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]

A, B, C, D = sp.symbols('A B C D')
GEN = [A, B, C, D]
def formalize(S, c0):
    out = {}
    for mono, v in S.items():
        nm = []
        for t in mono:
            if sp.simplify(t - 2) == 0: nm.append(A)
            elif sp.simplify(t - c0) == 0: nm.append(B)
            elif sp.simplify(t - (1 - c0)) == 0: nm.append(C)
            elif sp.simplify(t - (1 + c0)) == 0: nm.append(D)
            else:
                raise ValueError(f"unexpected letter {t}")
        out[tuple(nm)] = out.get(tuple(nm), 0) + v
    return {k: v for k, v in out.items() if v}

def sym_Omega_formal(c0):
    S = {}
    for w, cf in TERMS:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf))
    return formalize(S, c0)

# ---------- embeddings & relations ----------
def embeddings(c0):
    """all real embeddings of Q(c0) as numerical values of c0."""
    mpoly = sp.minpoly(c0, sp.Symbol('xx'))
    xx = sp.Symbol('xx')
    roots = sp.real_roots(sp.Poly(mpoly, xx))
    vals = []
    for r in roots:
        vals.append(mp.mpf(str(r.evalf(60))))
    return mpoly, vals

def letter_vals(c0val):
    return [mp.mpf(2), c0val, 1 - c0val, 1 + c0val]

def log_matrix(c0):
    mpoly, cvals = embeddings(c0)
    M = []
    for cv in cvals:
        row = [mp.log(abs(v)) for v in letter_vals(cv)]
        M.append(row)
    return mpoly, M   # d embeddings x 4 letters

from lll_rel import find_relations

def find_relations_exact(c0):
    """find integer relations among [2],[c0],[1-c0],[1+c0] via log-embedding LLL;
       verify each exactly in the field; return verified relation lattice (rows)."""
    mpoly, M = log_matrix(c0)
    # lattice: letters = 4 columns; relations are integer vectors e with sum_i e_i log|sigma_j g_i| = 0 for all j
    # use the standard trick: find_relations on the concatenated log values
    xs = []
    for j, row in enumerate(M):
        xs.extend(row)
    # We want vectors e in Z^4 with M e = 0. Build lattice directly:
    # rows: (e_i, N*(M e_i)) -> kernel vectors
    N = mp.mpf(10) ** 60
    Bmat = []
    for i in range(4):
        r = [0] * (4 + len(xs))
        r[i] = 1
        for j, x in enumerate(xs):
            r[4 + j] = int(mp.nint(xs[j] * N))
        Bmat.append(r)
    import flint
    R = flint.fmpz_mat(Bmat).lll().tolist()
    rels = []
    for r in R[:4]:
        e = r[:4]
        # check residual
        res = [sum(e[i] * M[j][i] for i in range(4)) for j in range(len(M))]
        if max(abs(x) for x in res) < mp.mpf('1e-50'):
            rels.append(e)
    # exact verification: prod g_i^e_i == +-1 in the field
    g = [sp.Integer(2), c0, 1 - c0, 1 + c0]
    vrels = []
    for e in rels:
        prod = sp.Integer(1)
        for gi, ei in zip(g, e):
            if ei >= 0:
                prod *= gi ** ei
            else:
                prod /= gi ** (-ei)
        prod = sp.simplify(prod)
        if prod in (1, -1):
            vrels.append((e, prod))
        else:
            print("   WARNING: LLL relation not exact:", e, "prod =", prod)
    return mpoly, vrels, M

# ---------- delta22 in the presentation ----------
def wedge_basis_relations(rels):
    """given relation rows (e in Z^4: sum e_i [g_i] = 0), produce the induced subspace
       of Lambda^2(F4) (relations among wedges) and of Lambda^2(Lambda^2 F4)."""
    # F4 basis: A,B,C,D (indices 0..3). Relation r: sum r_i g_i = 0.
    pass

def delta22_verdict(T, rels):
    """T: formal delta22 dict over {A,B,C,D}; rels: verified relation rows.
       Compute image of T in Lambda^2(Lambda^2(F4/R)) and whether it's zero,
       where R = span(rels). Returns (nonzero_in_quotient, dim_U, details)."""
    import sympy as sp
    R = sp.Matrix([list(r) for r in rels]) if rels else sp.zeros(1, 4)
    if rels:
        ns = R.nullspace()
    else:
        ns = [sp.Matrix(v) for v in sp.eye(4).tolist()]
    # basis of U = F4/R: columns of nullspace(R)
    basisU = [list(v) for v in ns]
    k = len(basisU)
    # express each generator g_i in basis U-coords: g_i = sum_j m_ij u_j  (projection along basis)
    # The map F4 -> U: since basisU is a basis of U embedded in F4, coordinates of g_i:
    # not unique unless we pick a section; but the wedge map factors through F4 -> U canonically:
    # take the quotient map pi: F4 -> R^k via coordinates in any complement.
    # Simpler: work directly in Lambda^2(F4) and quotient by relations-induced subspace.
    # Lambda^2(F4): basis e_ij (i<j). Relations: for each r in R and each generator g_l:
    #    (sum_i r_i g_i) ^ g_l = sum_i r_i (e_i ^ g_l) = 0  -> relation among wedges.
    # quotient Lambda^2(U) = Lambda^2(F4)/span{those}.
    import itertools
    pairs = list(itertools.combinations(range(4), 2))
    pidx = {p: q for q, p in enumerate(pairs)}
    rel_rows = []
    for r in rels:
        for l in range(4):
            row = [0] * 6
            for i in range(4):
                if i == l:
                    continue
                p = (min(i, l), max(i, l))
                sgn = 1 if i < l else -1
                row[pidx[p]] += r[i] * sgn
            rel_rows.append(row)
    # also the quadratic relations: (sum r_i e_i) ^ (sum s_j e_j) for r,s in R? those are in span already
    # Lambda^2(U) basis: quotient; then Lambda^2(Lambda^2 U): quotient of Lambda^2(6-dim) by
    # relations (wedge relations among the 6 pairs) induced by wedge-relations.
    Rw = sp.Matrix(rel_rows) if rel_rows else sp.zeros(1, 6)
    ns2 = Rw.nullspace()
    basisW = [list(v) for v in ns2]   # basis of Lambda^2(U) as combos of e_ij
    kW = len(basisW)
    # T in coordinates of e_ij-wedges (first-level wedge): T is sum over ((p),(q)) of
    # coeff * e_p ^ e_q (outer). Map T to Lambda^2(Lambda^2 U):
    # express each e_p in basisW coords (projection): solve for the induced map.
    # Build the linear map Lambda^2(F4) -> Lambda^2(U) = span(basisW): coordinates via
    # the quotient: take matrix Wm with columns basisW; coords of a vector v in quotient
    # = least solution; instead: quotient map to R^{6}/im(Rw): use projection onto a complement.
    # Practical: pick a basis of the quotient via rref of Rw: free columns.
    if kW == 0:
        return False, 0, {"k": k}
    # quotient coordinates: choose rows of a left inverse: compute row-echelon of Rw and
    # use the "free coordinate" functionals.
    # Simpler robust approach: build matrix Q whose rows = basisW^T (as combos of e_ij),
    # map v -> (coords w.r.t. basisW via solving Wm * c = v + rel-adjustment):
    # Since basisW spans the quotient, and any v in R^6 maps: v mod im(Rw).
    # Coordinates: choose a set of columns of identity forming a complement basis.
    RwT = Rw.T
    # find column basis of RwT (image): pivot columns
    M = sp.Matrix.hstack(RwT, sp.eye(6))
    rrefM, piv = M.rref()
    # complement coordinates: coordinates along non-pivot columns of the identity part? Not directly.
    # Use: quotient map given by [0 .. I_kW] after column operations:
    # Compute GL: invertible 6x6 P with P*Rw in rref; quotient = last kW rows of P.
    P, pivots = (sp.Matrix.vstack(Rw, sp.eye(6))).rref()  # not the way
    # do it concretely: solve for coordinates of T's first-level wedges inside span(basisW)
    # T's outer wedge: sum over (p,q) t_{pq} e_p ^ e_q.
    # Compute image: for each e_p, its class [e_p] in quotient; the outer products [e_p]^[e_q].
    # The induced relations on the quotient's wedges: [e_p]^[e_q] with relations from
    # basis of Lambda^2(U) (dim kW). The quotient Lambda^2(U) has its own wedge pairing
    # into Lambda^2(Lambda^2 U), dim kW choose 2. Map e_p (6 pairs) -> quotient (kW basis):
    # coordinates: pick a section: matrix S (6 x kW) = basisW as columns; for each e_p,
    # write e_p = S c_p + rw (rw in im(Rw)); coords c_p by least squares exact solve.
    Sm = sp.Matrix.hstack(*[sp.Matrix(v) for v in basisW]) if basisW else sp.zeros(6, 1)
    # Solve Sm c_p = e_p + (element of im Rw): i.e., solve [Sm | Rw] [c; s] = e_p
    aug = sp.Matrix.hstack(Sm, Rw)
    def coords(v):
        sol = aug.gauss_jordan_solve(sp.Matrix(v)) if aug.cols else None
        try:
            s = aug.gauss_jordan_solve(sp.Matrix(v))
            return [s[i] for i in range(Sm.cols)]
        except Exception:
            return None
    ecoords = [coords(sp.Matrix([1 if j == q else 0 for j in range(6)])) for q in range(6)]
    if any(c is None for c in ecoords):
        return None, k, {"error": "coord failure"}
    # T's value in Lambda^2(Lambda^2 U): sum_{p,q} t_{pq} (sum_a c_{p,a} u_a) ^ (sum_b c_{q,b} u_b)
    # outer wedge antisymmetric in the two blocks:
    Tsum = {}
    for (p, q), v in T.items():
        # p, q are pairs of generators ((g1,g2),(g3,g4)) as wedge keys; map to 6-dim pair indices
        pass
    return None, k, {"kW": kW}

print("computing formal symbols at all points...")
POINTS = [(f"sin(pi/{n})", sp.sin(sp.pi/n)) for n in range(3, 13)] + [("cos(pi/12)", sp.cos(sp.pi/12))]
formal = {}
for name, c0 in POINTS:
    formal[name] = sym_Omega_formal(c0)
    print(name, "symbol terms:", len(formal[name]), " delta22 terms:", len(delta22(formal[name])))
import pickle
pickle.dump({k: (str(v), str(delta22(v))) for k, v in formal.items()},
            open("<project>/research/08_theory/COR6X_independent_corner_symbols/formal_symbols.pkl", "wb"))
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_08.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_08.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout[-3000:]); print("ERR:", (r.stderr or "")[-2000:])
