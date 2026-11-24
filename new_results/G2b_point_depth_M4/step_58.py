
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
from cor6_sym import sshuffle, rref
from itertools import product
C0 = Fr(1,2)
import g2_fit as gf
from step_10 import load_terms

# shuffle subspace D^3 over 2 letters (exponent-vector basis letters e0=(1,0), e1=(0,1))
e0, e1 = (Fr(1), Fr(0)), (Fr(0), Fr(1))
Drows = []
B3 = list(product([e0, e1], repeat=3))
for a in (e0, e1):
    for bc in product([e0, e1], repeat=2):
        sh = sshuffle({(a,): Fr(1)}, {bc: Fr(1)})
        Drows.append([sh.get(w, Fr(0)) for w in B3])
Dmat, piv = rref([r[:] for r in Drows])
print("dim D^3:", len(Dmat))

def in_d3(S):
    v = [S.get(w, Fr(0)) for w in B3]
    M2, p2 = rref([r[:] for r in Drows] + [v])
    return len(M2) == len(Dmat)

words5 = sorted({w for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']) if len(w) == 5})
ok, bad = 0, 0
for w in words5:
    D5w = s05.delta_minus_5sym(s05.S_J(w, C0, 0))
    Dw = s05.delta_word_sym(w, 0, C0)
    # delta_word + delta^-(sym) should be a shuffle product (vanish in w3 mod products)
    diff = dict(Dw)
    for k, v in D5w.items():
        diff[k] = diff.get(k, Fr(0)) + v
    diff = {k: v for k, v in diff.items() if v}
    if in_d3(diff):
        ok += 1
    else:
        bad += 1
        if bad <= 2:
            print("NOT in D^3:", w, diff)
print(f"delta_word + delta^-(sym) in D^3 for {ok}/{ok+bad} words")
