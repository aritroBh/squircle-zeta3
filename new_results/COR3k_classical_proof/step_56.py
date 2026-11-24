import json, re, sympy as sp
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
c = sp.symbols('c', positive=True)
q = (1 - c**2)/(2*c)
R2 = json.load(open(G1 + 'hpl_fit_k2_(0,0,0).json'))

# ---- typed B2 (result.md): M2 = z2/2 - Li2(c^2)/2 - z2*q*artanh(c) + S(c) + S(-c)
#      S(c) = (1-c)[H001 + H01-1 - H1-11] + q H101 - ((1-c)^2/(2c)) H10-1
# Expand S(c)+S(-c) via the sign rule H_w(-c) = (-1)^(#nonzero letters) H_{-w}(c):
def nnz(w): return sum(1 for a in w if a != 0)
typed_hpl = {}
def addh(w, cf):
    typed_hpl[tuple(w)] = sp.simplify(typed_hpl.get(tuple(w), 0) + cf)
Sc_terms = [((0,0,1), (1-c)), ((0,1,-1), (1-c)), ((1,-1,1), -(1-c)),
            ((1,0,1), q), ((1,0,-1), -(1-c)**2/(2*c))]
for w, cf in Sc_terms:                       # S(c)
    addh(w, cf)
for w, cf in Sc_terms:                       # S(-c): c -> -c, H_w(-c) -> (-1)^nnz H_{-w}(c)
    addh(tuple(-a for a in w), (-1)**nnz(w) * sp.expand(cf.subs(c, -c)))
addh((0, 1), -1)      # -Li2(c^2)/2 = -Li2(c) - Li2(-c)   (duplication, proved below)
addh((0, -1), 1)

json_hpl = {}
for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = tuple(int(z) for z in mm.group(2).split(','))
    json_hpl[w] = sp.simplify(json_hpl.get(w, 0) + sp.Rational(cstr)*c**e)

print("=== typed B2  vs  hpl_fit_k2_(0,0,0).json  (coefficients in Q[c,1/c]) ===")
allok = True
for w in sorted(set(typed_hpl) | set(json_hpl)):
    d = sp.simplify(typed_hpl.get(w, 0) - json_hpl.get(w, 0))
    allok &= (d == 0)
    print(f"  H{list(w)}: typed {sp.factor(typed_hpl.get(w,0))}  vs  json {sp.factor(json_hpl.get(w,0))}  diff {d}")
print("json == typed B2 (S(c)+S(-c) with sign rule):", allok)

# ---- duplication Li2(x) + Li2(-x) = Li2(x^2)/2 (used in B2 -> B2') ----
x2 = sp.symbols('x2')
Gdup = sp.polylog(2, x) + sp.polylog(2, -x) - sp.polylog(2, x**2)/2
dG = sp.expand(dpoly(Gdup))
okd, _ = check_weight1_zero(dG)
print("\nduplication: derivative log-coeffs all 0:", okd, "; value at 0 is 0 (Li2(0)=0) ->",
      "PROVED" if okd else "FAILED")
