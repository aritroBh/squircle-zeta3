import sympy as sp
from cor3_common import *

# duplication Li2(x) + Li2(-x) = Li2(x^2)/2:
# d/dx by hand: polylog(1, x^2) -> -log(1-x^2) = -(log(1-x) + log(1+x))  (both factors > 0 on (0,1))
Gdup = sp.polylog(2, x) + sp.polylog(2, -x) - sp.polylog(2, x**2)/2
dG = sp.diff(Gdup, x)
dG = dG.replace(sp.polylog(1, x**2), -(l1m + l1p))   # hand factorization, valid on (0,1)
dG = replace_polylog1(dG)
okd, terms = check_weight1_zero(sp.expand(dG))
c0, _ = const_term(sp.polylog(2, x) + sp.polylog(2, -x), 2)   # value at 0 = 0 (Li2(0)=0)
print("duplication: d/dx[LHS-RHS] log-coeffs:", terms, "all zero:", okd,
      "; value at 0+:", c0, "->", "PROVED" if okd and c0 == 0 else "FAILED")

# ---- sign rule H_w(-x) = (-1)^nnz H_{-w}(x): numeric spot check (proof is one line, see report)
import mpmath as mp
mp.mp.dps = 50
def f_at(a, t):
    return 1/t if a == 0 else (1/(1-t) if a == 1 else 1/(1+t))
def Hw(w, X):
    if len(w) == 1:
        a = w[0]
        return mp.log(X) if a == 0 else (-mp.log(1-X) if a == 1 else mp.log(1+X))
    return mp.quad(lambda t: f_at(w[0], t)*Hw(w[1:], t), [0, X])
for w in [(0,0,1), (0,1,-1), (1,-1,1), (1,0,-1)]:
    for X in [mp.mpf('0.3'), mp.mpf('0.7')]:
        lhs = Hw(w, -X); rhs = (-1)**sum(1 for a in w if a != 0)*Hw(tuple(-a for a in w), X)
        print(f"  sign rule {w} x={X}: diff {mp.nstr(abs(lhs-rhs), 2)}")
