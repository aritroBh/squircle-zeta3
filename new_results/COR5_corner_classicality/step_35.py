import sys, json
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

# 21 S-unit solutions
sols = json.load(open("sunit23_solutions.json"))
def vec23(fac):
    s, a, b = fac
    return (Fr(a), Fr(b))   # sign ignored: [-1]=0 in letter space
li4 = []
for xs, fx, fy in sols:
    x = Fr(xs)
    rv, om = vec23(fx), vec23(fy)
    S = li_sym(4, rv, om)
    v = Q.coords(S)
    li4.append((x, v))
print("projected Li4 symbols (quotient coords):")
for x, v in li4:
    print(f"  x = {str(x):>6}: {v}")

# rank of all 21
M = [list(v) for _, v in li4]
Dm, piv = rref(M)
print("\nrank of all Li4 projected symbols:", len(Dm))
print("target pi(F(1/2)):", vF)
aug = rref([list(v) for _, v in li4] + [list(vF)])
print("rank with target adjoined:", len(aug), "=> in span:", len(aug) == len(Dm))
