import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *

Q = W4Quotient(2)
e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))

# HPL symbol of H(0,0,-1,1) at c=1/2: (-[1-c])[1+c][c][c]
HPL = {( (Fr(1),Fr(0)), (Fr(-1),Fr(1)), (Fr(-1),Fr(0)), (Fr(-1),Fr(0)) ): Fr(1)}
print("class(HPL symbol):", Q.coords(HPL))

# full identity RHS symbol:
# 6*Li4(-2): S(Li4(-2)) = -(e2,e1,e1,e1) -> 6*(-1)*(e2,e1,e1,e1)
# -2*Li4(1/2): S(Li4(1/2)) = -(-e1,-e1,-e1,-e1) = +(e1,e1,e1,e1) -> -2*(e1,e1,e1,e1)
S = {}
def add(S0, q):
    for w, c in S0.items():
        S[w] = S.get(w, 0) + q*c
add({((e2,e1,e1,e1)): Fr(-6)}, Fr(1))
add({((e1,e1,e1,e1)): Fr(-2)}, Fr(1))
# -Li3(-1/2)*log2: S(Li3(-1/2)) = -((-e1+e2),(-e1),(-e1)) ; shuffle with (e1), coeff -1
S_Li3mh = {(ladd(lneg(e1),e2), lneg(e1), lneg(e1)): Fr(-1)}
add(sshuffle(S_Li3mh, {(e1,): Fr(1)}), Fr(-1))
# -(3/4)zeta3*log2, (5/8)pi2*log2^2, (137/1440)pi4: symbol 0 (zeta3, pi2, pi4 invisible)
# (1/24)log2^4: (e1)^4 shuffle -> in D^4 anyway
print("class(identity RHS full symbol):", Q.coords(S))
print("difference:", tuple(a-b for a, b in zip(Q.coords(HPL), Q.coords(S))))
