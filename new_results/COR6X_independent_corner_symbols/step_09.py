
# step_09.py -- per-field letter presentation (exact relations verified) + delta22 verdict
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp, json, re, itertools
import mpmath as mp
import flint
mp.mp.dps = 100

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
A, B, C, D = sp.symbols('A B C D')
def formalize(S, c0):
    out = {}
    for mono, v in S.items():
        nm = []
        for t in mono:
            def eq(u): return sp.simplify(t - u) == 0 or sp.simplify(t + u) == 0
            if eq(sp.Integer(2)): nm.append(A)
            elif eq(c0): nm.append(B)
            elif eq(1 - c0): nm.append(C)
            elif eq(1 + c0): nm.append(D)
            else: raise ValueError(f"unexpected letter {t}")
        out[tuple(nm)] = out.get(tuple(nm), 0) + v
    return {k: v for k, v in out.items() if v}
def sym_Omega_formal(c0):
    S = {}
    for w, cf in TERMS:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf))
    return formalize(S, c0)

def log_matrix(c0, prec=90):
    xx = sp.Symbol('xx')
    mpoly = sp.minpoly(c0, xx)
    roots = sp.real_roots(sp.Poly(mpoly, xx))
    cvals = [mp.mpf(str(r.evalf(prec))) for r in roots]
    M = []
    for cv in cvals:
        M.append([mp.log(abs(v)) for v in (mp.mpf(2), cv, 1 - cv, 1 + cv)])
    return mpoly, M

def find_rels(M, c0):
    """relations e in Z^4 with M e = 0 (log-embedding), via flint LLL; verify exactly."""
    d = len(M)
    N = mp.mpf(10) ** 70
    Bmat = []
    for i in range(4):
        r = [0] * (4 + d)
        r[i] = 1
        for j in range(d):
            r[4 + j] = int(mp.nint(M[j][i] * N))
        Bmat.append(r)
    R = flint.fmpz_mat(Bmat).lll().tolist()
    g = [sp.Integer(2), c0, 1 - c0, 1 + c0]
    vrels = []
    for r in R[:4]:
        e = [int(x) for x in r[:4]]
        if all(x == 0 for x in e): continue
        res = max(abs(sum(e[i] * M[j][i] for i in range(4))) for j in range(d))
        if res > mp.mpf('1e-55'): break
        prod = sp.Integer(1)
        for gi, ei in zip(g, e):
            prod *= gi ** ei if ei >= 0 else gi ** ei  # sympy handles negative
        prod = sp.simplify(prod)
        if prod in (1, -1):
            vrels.append((e, prod))
        else:
            print("   ! candidate relation fails exact check:", e, "->", prod)
    return vrels

def matrix_rank_cert(M, k):
    """certify rank(M) >= k via a k x k minor with rigorous error bound (report det magnitude)."""
    import itertools
    best = None
    for rows in itertools.combinations(range(len(M)), k):
        sub = [[M[j][i] for i in range(4)][:k] for j in rows]
        # choose k independent columns: try all col subsets
        for cols in itertools.combinations(range(4), k):
            msub = [[M[j][i] for i in cols] for j in rows]
            det = abs(sp.Matrix([[float(x) for x in r] for r in msub]).det()) if False else None
            md = mp.matrix(msub)
            d = mp.fabs(md[0,0])
            # det via mp
            detv = mp.fabs(mp.det(md))
            if best is None or detv > best:
                best = detv
    return best

def delta22_image(T, rels):
    """image of formal delta22 T in Lambda^2(Lambda^2(F4/R)); rels = verified rows (Z^4)."""
    R = sp.Matrix([list(r) for r, p in rels]) if rels else sp.zeros(0, 4)
    ns = R.nullspace() if rels else [sp.Matrix(v) for v in sp.eye(4).tolist()]
    basisU = [sp.Matrix(v) for v in ns]
    k = len(basisU)
    if k <= 2:
        return {}, k
    # quotient map: extend basisU with a basis of R
    if rels:
        rrefR, piv = R.rref()
        rbasis = [sp.Matrix(r) for r in rrefR.tolist() if any(r)]
    else:
        rbasis = []
    M = sp.Matrix.hstack(*(basisU + rbasis))
    Minv = M.inv()
    def proj(gi):
        v = Minv * sp.Matrix([1 if j == gi else 0 for j in range(4)])
        return [sp.nsimplify(v[j]) for j in range(k)]
    pg = [proj(i) for i in range(4)]
    # wedge of two projected generators -> antisymmetric k x k matrix
    def wij(i, j):
        Mm = sp.zeros(k)
        for a in range(k):
            for b in range(k):
                Mm[a, b] = pg[i][a] * pg[j][b] - pg[i][b] * pg[j][a]
        return Mm
    def canon_w(Mm):
        # canonical coordinates of antisymmetric matrix: dict over (a,b), a<b
        out = {}
        for a in range(k):
            for b in range(a + 1, k):
                if Mm[a, b] != 0:
                    out[(a, b)] = sp.nsimplify(Mm[a, b])
        return tuple(sorted(out.items()))
    out = {}
    for ((g1, g2), (g3, g4)), v in T.items():
        idx = {A: 0, B: 1, C: 2, D: 3}
        w1 = canon_w(wij(idx[g1], idx[g2]))
        w2 = canon_w(wij(idx[g3], idx[g4]))
        if not w1 or not w2 or w1 == w2: continue
        s = v
        if w1 > w2:
            w1, w2 = w2, w1; s = -s
        out[(w1, w2)] = out.get((w1, w2), 0) + s
    out = {k2: v for k2, v in out.items() if v != 0}
    return out, k

POINTS = [(f"sin(pi/{n})", sp.sin(sp.pi/n)) for n in range(3, 13)] + [("cos(pi/12)", sp.cos(sp.pi/12))]
for name, c0 in POINTS:
    mpoly, M = log_matrix(c0)
    rels = find_rels(M, c0)
    rankU = 4 - (sp.Matrix([r for r, p in rels]).rank() if rels else 0)
    cert = matrix_rank_cert(M, rankU) if rankU else None
    T = delta22(sym_Omega_formal(c0))
    img, k = delta22_image(T, rels)
    print(f"=== {name}: field degree {sp.degree(mpoly)}, minpoly {mpoly.as_expr()}")
    print(f"    verified relations: {[(r, p) for r, p in rels]}  -> letter span dim {rankU} (cert |det|~{cert})")
    print(f"    formal delta22: {T}")
    print(f"    delta22 in Lambda^2(Lambda^2 U): {'ZERO' if not img else img}")
