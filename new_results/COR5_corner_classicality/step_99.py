import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
from itertools import product as iproduct

# weight-3 quotient Q^3 = T^3/D^3 for V = Q^2
e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
letters = [e1, e2]
W3 = list(iproduct(letters, repeat=3))
W1 = list(iproduct(letters, repeat=1))
W2 = list(iproduct(letters, repeat=2))
D = []
for u in W1:
    for v in W2: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), 2))
for u in W2:
    for v in W1: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), 2))
def vec_of(S): return tuple(S.get(w, Fr(0)) for w in W3)
Dm, piv = rref([vec_of(s) for s in D])
free = [j for j in range(8) if j not in piv]
def q3(S):
    v = list(vec_of(expand_sym(S, 2)))
    for i, pj in enumerate(piv):
        if v[pj]:
            q = v[pj]
            v = [a - q*b for a, b in zip(v, Dm[i])]
    assert all(v[j] == 0 for j in range(8) if j not in free)
    return tuple(v[j] for j in free)
print("free words:", [ "".join('2' if l==e1 else '3' for l in W3[j]) for j in free])
S_Li3_m2 = li_sym(3, e1, e2)   # S(Li3(-2)): rv=[-2]=e1, om=[3]=e2
print("S(Li3(-2)):", S_Li3_m2, " class:", q3(S_Li3_m2))
S_Li3_half = li_sym(3, lneg(e1), lneg(e1))  # S(Li3(1/2)): rv=[1/2]=-e1, om=[1/2]=-e1
print("S(Li3(1/2)) class:", q3(S_Li3_half))
S_Li3_mhalf = li_sym(3, lneg(e1), ladd(lneg(e1), e2))  # S(Li3(-1/2)): rv=[1/2]=-e1, om=[3/2]=-e1+e2
print("S(Li3(-1/2)) class:", q3(S_Li3_mhalf))
# identity's weight-3 combination: 6Li3(-2) + 2Li3(1/2) - Li3(-1/2) - (3/4)zeta3
comb = sadd(sscale(6, S_Li3_m2), sscale(2, S_Li3_half), sscale(-1, S_Li3_mhalf))
print("identity weight-3 class:", q3(comb))
# HPL's delta_31 Lie part: -Lie(S(Li3(-2))) per my hand computation
print("HPL Lie part class = -S(Li3(-2)) class:", tuple(-x for x in q3(S_Li3_m2)))
