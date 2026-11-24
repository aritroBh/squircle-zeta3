import sys, json, itertools
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
import sympy as sp

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
vecs_all = [v for _, v in li4]; xs_all = [x for x, _ in li4]

def solve_span(vecs, target):
    A = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in v] for v in vecs]).T
    b = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in target])
    if A.rank() != A.row_join(b).rank(): return None
    sol = A.gauss_jordan_solve(b)
    par = sol[0]
    for p in sol[1]: par = par.subs(p, 0)
    return [Fr(par[i].p, par[i].q) for i in range(len(vecs))]

def height(cs):
    return max(max(abs(c.numerator), c.denominator) for c in cs)

# all 3-subset solutions, sorted by height
sols3 = []
for idxs in itertools.combinations(range(len(li4)), 3):
    cs = solve_span([vecs_all[i] for i in idxs], vF)
    if cs is not None:
        sols3.append((height(cs), idxs, cs))
sols3.sort()
print("10 lowest-height 3-argument solutions:")
for h, idxs, cs in sols3[:10]:
    print(f"  height {h:>4}: args {[str(xs_all[i]) for i in idxs]}  coeffs {[str(c) for c in cs]}")
print()
print("solutions with all args in (0,1):")
for h, idxs, cs in sols3:
    if all(0 < xs_all[i] < 1 for i in idxs):
        print(f"  height {h:>4}: args {[str(xs_all[i]) for i in idxs]}  coeffs {[str(c) for c in cs]}")
