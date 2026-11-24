import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
from itertools import product as iproduct

# Eulerian idempotent e^{(1)} on words over basis letters of V=Q^2:
# e(w) = sum_{k>=1} (-1)^{k-1}/k * sum_{(u1..uk) |- w, ui nonempty} u1 ш ... ш uk
def eulerian(w):
    out = {}
    n = len(w)
    # compositions of the word into k consecutive nonempty blocks
    def comps(w):
        # yield all tuples of consecutive nonempty subwords
        if not w:
            yield ()
            return
        for i in range(1, len(w)+1):
            for rest in comps(w[i:]):
                yield (w[:i],) + rest
    for comp in comps(w):
        k = len(comp)
        coef = Fr((-1)**(k-1), k)
        S = {comp[0]: Fr(1)}
        for blk in comp[1:]:
            S = sshuffle(S, {blk: Fr(1)})
        for u, c in S.items():
            out[u] = out.get(u, 0) + coef * c
    return {u: c for u, c in out.items() if c}

e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
# test: kills shuffles
t1 = eulerian((e1, e1, e2, e1))
Ssh = sshuffle({(e1,): Fr(1)}, {(e1, e2, e1): Fr(1)})
tot = {}
for w, c in Ssh.items():
    for u, cc in eulerian(w).items():
        tot[u] = tot.get(u, 0) + c*cc
tot = {u: c for u, c in tot.items() if c}
print("eulerian kills shuffle (e1)ш(e1 e2 e1):", tot == {})
# test: Li4 symbol survives
Sli4 = li_sym(4, e1, e2)   # -[1-x][x]^3 with rv=e1, om=e2 -> word (-e2, e1,e1,e1)
tot2 = {}
for w, c in Sli4.items():
    # expand -e2 = -1 * e2
    for u, cc in eulerian(w).items():
        tot2[u] = tot2.get(u, 0) + c*cc
tot2 = {u: c for u, c in tot2.items() if c}
print("eulerian(S(Li4)) nonzero:", tot2 != {})
# compare with quotient: eulerian coords on the free basis words of Q
Q = W4Quotient(2)
freewords = [Q.W4[j] for j in Q.free]
def evec(S):
    tot = {}
    for w, cf in expand_sym(S, 2).items():
        for u, cc in eulerian(w).items():
            tot[u] = tot.get(u, 0) + cf*cc
    tot = {u: c for u, c in tot.items() if c}
    return tuple(tot.get(w, Fr(0)) for w in freewords)
print("eulerian-coords(S(Li4) rv=e1,om=e2):", evec(Sli4), " quotient:", Q.coords(Sli4))
