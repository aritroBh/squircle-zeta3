import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
import random

# --- validate the Dynkin projection ---
# V = Q{e1, e2}; letters as 2-dim vectors
e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
letters = [e1, e2]

# 1) image span of proj_lie on all 16 words is 3-dimensional (= dim Lie_4 for 2 letters)
from itertools import product as iproduct
def words_of_len(n, lets):
    return list(iproduct(lets, repeat=n))

W4 = words_of_len(4, letters)
def vec_of(L):
    # L: dict word->Fr over the 16 words -> tuple over fixed word list
    return tuple(L.get(w, Fr(0)) for w in W4)

def rank(rows):
    rows = [list(r) for r in rows]
    R = 0; piv = []
    for j in range(len(rows[0])):
        p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
        if p is None: continue
        rows[R], rows[p] = rows[p], rows[R]
        q = rows[R][j]
        rows[R] = [x/q for x in rows[R]]
        for i in range(len(rows)):
            if i != R and rows[i][j]:
                rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
        piv.append(j); R += 1
        if R == len(rows): break
    return R, piv, rows

imgs = [vec_of(dynkin_left(w)) for w in W4]
R, piv, RR = rank(imgs)
print("rank of Dynkin images on 16 words:", R, "(expect 3)")

# 2) Ree/Dynkin-Specht-Wever: proj kills shuffle products  (random checks + exhaustive small)
random.seed(1)
ok = True
# exhaustive: all splits of words into (u,v) with |u|+|v|=4, letters e1/e2, random coef symbols
tested = 0
for nu in (1, 2, 3):
    nv = 4 - nu
    for u in words_of_len(nu, letters):
        for v in words_of_len(nv, letters):
            S = sshuffle({u: Fr(1)}, {v: Fr(1)})
            P = proj_lie(S)
            if P:
                ok = False
                print("FAIL:", u, v, P)
                break
            tested += 1
print(f"shuffle products killed by proj_lie: {tested} exhaustive checks, all zero: {ok}")

# 3) Lie elements are preserved up to factor 4: check dynkin_left(P) == 4*pi... on Lie basis
# Lyndon basis elements: P[1112] = [e1,[e1,[e1,e2]]], etc.
def lyndon_bracket(*elts):
    # build from nested bracket spec
    pass
def Lbr(A, B):
    out = {}
    for u, cu in A.items():
        for v, cv in B.items():
            for w, c in lieb_words(u, v).items():
                out[w] = out.get(w, 0) + cu*cv*c
    return {w: c for w, c in out.items() if c}
E1 = {(e1,): Fr(1)}; E2 = {(e2,): Fr(1)}
L1 = Lbr(E1, Lbr(E1, Lbr(E1, E2)))          # [e1,[e1,[e1,e2]]]
L2 = Lbr(E1, Lbr(Lbr(E1, E2), E2))          # [e1,[[e1,e2],e2]]
L3 = Lbr(Lbr(Lbr(E1, E2), E2), E2)          # [[[e1,e2],e2],e2]
for name, L in [("L1", L1), ("L2", L2), ("L3", L3)]:
    # proj of the weight-4 symbol L (which is a Lie element): pi(L) should equal L
    P = proj_lie(L)
    print(name, "proj(L)==L:", P == L, " #terms:", len(L))
