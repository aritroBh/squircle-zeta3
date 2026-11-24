import sys
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor5_sym import *
import cor4_sym as c4

e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
# HPL symbol of H(0,0,-1,1) at c=1/2:
HPL = {((Fr(1),Fr(0)), (Fr(-1),Fr(1)), (Fr(-1),Fr(0)), (Fr(-1),Fr(0))): Fr(1)}
d31_HPL = c4.delta31(HPL, 2)
print("delta_31(HPL symbol):")
for k, v in sorted(d31_HPL.items()):
    print("   ", k, v)

# identity RHS weight-4 symbols:
# 6*Li4(-2): S = 6*(-[3][2]^3) = -6(e2,e1,e1,e1)
# -2*Li4(1/2): S = -2*(-(-e1)^4) = +2(e1,e1,e1,e1)
# -Li3(-1/2)*log2: S = -{(-[3/2][1/2]^2) ш (e1)}
# -(3/4)zeta3*log2: 0
# products: 0 in delta31? pi2*log2^2: S=0; log2^4: shuffle; pi4: 0
S_id = {}
for w, c in sscale(-6, {( (e2,e1,e1,e1) ): Fr(1)}).items(): S_id[w] = S_id.get(w,0)+c
for w, c in sscale(2, {((e1,e1,e1,e1)): Fr(1)}).items(): S_id[w] = S_id.get(w,0)+c
# -Li3(-1/2)*log2: S(Li3(-1/2)) = -[3/2][1/2]^2 = -((-e1+e2),(-e1),(-e1)) word
S_Li3mh = {(ladd(lneg(e1),e2), lneg(e1), lneg(e1)): Fr(-1)}
S_prod = sscale(-1, sshuffle(S_Li3mh, {(e1,): Fr(1)}))
for w, c in S_prod.items(): S_id[w] = S_id.get(w,0)+c
S_id = {w: c for w, c in S_id.items() if c}
d31_id = c4.delta31(S_id, 2)
print("delta_31(identity RHS):")
for k, v in sorted(d31_id.items()):
    print("   ", k, v)
