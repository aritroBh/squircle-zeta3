import sys, json
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
M = [list(v) for _, v in li4] + [list(vF)]
print("vF == (5/2,1/2,0):", vF == (Fr(5,2), Fr(1,2), Fr(0)))
As = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in M])
print("sympy rank of 22x3:", As.rank())
# now my rref, with prints
rows = [list(r) for r in M]; R = 0; piv = []
for j in range(len(rows[0])):
    p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
    print(f"col {j}: R={R} pivot={p}")
    if p is None: continue
    rows[R], rows[p] = rows[p], rows[R]
    q = rows[R][j]; rows[R] = [x/q for x in rows[R]]
    for i in range(len(rows)):
        if i != R and rows[i][j]:
            rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
    piv.append(j); R += 1
    if R == len(rows): break
print("R =", R)
# print the li4 vectors again to compare with the earlier debug cell
for x, v in li4:
    print(x, v)
