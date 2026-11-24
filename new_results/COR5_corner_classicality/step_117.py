import sympy as sp
from sympy import sin, pi
from sympy.polys.numberfields import minimal_polynomial
import mpmath as mp
from flint import fmpz_mat, fmpz

def log_vectors(elements, prec=200):
    """elements: list of sympy algebraics. Return list of log|sigma(el)| vectors over all conjugates."""
    mp.mp.dps = prec
    # use the field of the first element: conjugates via minpoly roots of the generator
    out = []
    for el in elements:
        x = sp.Symbol('x')
        P = minimal_polynomial(el, x)
        coeffs = [mp.mpf(sp.N(c, prec)) for c in sp.Poly(P).all_coeffs()]
        roots = mp.polyroots(coeffs, maxsteps=300)
        out.append([mp.log(abs(r)) for r in roots])
    return out

def find_relations(vecs, prec=200):
    """find integer relations among the log-vectors: rows (e_i | round(C * components))."""
    mp.mp.dps = prec
    n = len(vecs)
    dim = len(vecs[0])
    C = mp.mpf(10) ** (prec - 40)
    M = fmpz_mat(n, n + dim)
    for i, v in enumerate(vecs):
        M[i, i] = fmpz(1)
        for j, x in enumerate(v):
            M[i, n + j] = fmpz(int(mp.nint(C * x)))
    R = M.lll(delta=0.99, eta=0.51)
    rows = []
    for r in range(n):
        row = [int(R[r, c]) for c in range(n)]
        defect = [R[r, n + j] for j in range(dim)]
        rows.append((row, defect))
    return rows

# test on n=3: c = sqrt(3)/2, 1-c, 1+c
n = 3
c = sp.sqrt(3)/2
els = [c, 1 - c, 1 + c]
vecs = log_vectors(els)
rels = find_relations(vecs)
for row, defect in rels:
    a = row
    val = els[0]**a[0] * els[1]**a[1] * els[2]**a[2]
    val_s = sp.simplify(val)
    print(f"relation {a}: defect norm ~ {max(abs(int(d)) for d in defect)}, value = {val_s}")
