import sys, json, itertools
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *

base = "<project>/research"
cor4 = base + "/08_theory/COR4_M3alpha_classical"
Fterms = parse_F_symbol(cor4 + "/F_symbol.json")
letter_map = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(1))]
F_half = {}
for e, word, cf in Fterms:
    cval = Fr(2) ** (-e)
    nw = tuple(lincomb(l, letter_map) for l in word)
    F_half[nw] = F_half.get(nw, Fr(0)) + cf * cval
F_half = {w: c for w, c in F_half.items() if c}
Q = W4Quotient(2)
vF = Q.coords(F_half)
sols = json.load(open("sunit23_solutions.json"))
def vec23(fac):
    s, a, b = fac
    return (Fr(a), Fr(b))
li4 = [(Fr(xs), Q.coords(li_sym(4, vec23(fx), vec23(fy)))) for xs, fx, fy in sols]

def solve_span(vecs, target):
    """find coefficients c_i (Fractions) with sum c_i vecs[i] = target, or None."""
    import sympy as sp
    A = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in v] for v in vecs]).T
    b = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in target])
    if A.rank() != A.row_join(b).rank():
        return None
    par, params = A.gauss_jordan_solve(b)
    if params.rows:
        # set free params to 0
        sol = par
        for p in params:
            sol = sol.subs(p, 0)
    else:
        sol = par
    out = []
    for i in range(len(vecs)):
        r = sol[i]
        out.append(Fr(r.p, r.q))
    return out

# minimal combinations: try subsets of increasing size
vecs_all = [v for _, v in li4]
xs_all = [x for x, _ in li4]
def in_span(vecs, target):
    import sympy as sp
    A = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in v] for v in vecs]).T
    b = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in target])
    return A.rank() == A.row_join(b).rank()

found = None
for k in range(1, 6):
    goods = []
    for idxs in itertools.combinations(range(len(li4)), k):
        if in_span([vecs_all[i] for i in idxs], vF):
            goods.append(idxs)
    if goods:
        found = goods
        print(f"minimal number of Li4 arguments: {k}; number of minimal solutions: {len(goods)}")
        break
for idxs in found[:12]:
    cs = solve_span([vecs_all[i] for i in idxs], vF)
    print("  args:", [str(xs_all[i]) for i in idxs], " coeffs:", [str(c) for c in cs])
