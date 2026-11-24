import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
from itertools import product as iproduct

e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
letters = [e1, e2]
def words_of_len(n, lets): return list(iproduct(lets, repeat=n))
W4 = words_of_len(4, letters)
def vec_of(L): return tuple(L.get(w, Fr(0)) for w in W4)
def rank(rows):
    rows = [list(r) for r in rows]
    R = 0
    for j in range(len(rows[0])):
        p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
        if p is None: continue
        rows[R], rows[p] = rows[p], rows[R]
        q = rows[R][j]
        rows[R] = [x/q for x in rows[R]]
        for i in range(len(rows)):
            if i != R and rows[i][j]:
                rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
        R += 1
        if R == len(rows): break
    return R

imgs = [vec_of(dynkin_left(w)) for w in W4]
print("rank of Dynkin images on 16 words:", rank(imgs), "(expect 3)")

ok = True; tested = 0
for nu in (1, 2, 3):
    nv = 4 - nu
    for u in words_of_len(nu, letters):
        for v in words_of_len(nv, letters):
            S = sshuffle({u: Fr(1)}, {v: Fr(1)})
            if proj_lie(S):
                ok = False; print("FAIL:", u, v)
            tested += 1
print(f"Ree: shuffle products killed: {tested} exhaustive checks all zero: {ok}")

def Lbr(A, B):
    out = {}
    for u, cu in A.items():
        for v, cv in B.items():
            for w, c in lieb_words(u, v).items():
                out[w] = out.get(w, 0) + cu*cv*c
    return {w: c for w, c in out.items() if c}
E1 = {(e1,): Fr(1)}; E2 = {(e2,): Fr(1)}
L1 = Lbr(E1, Lbr(E1, Lbr(E1, E2)))
L2 = Lbr(E1, Lbr(Lbr(E1, E2), E2))
L3 = Lbr(Lbr(Lbr(E1, E2), E2), E2)
for name, L in [("L1", L1), ("L2", L2), ("L3", L3)]:
    print(name, "proj(L)==L:", proj_lie(L) == L, " #terms:", len(L))
# rank of the Lyndon basis Lie elements
print("rank of {L1,L2,L3}:", rank([vec_of(L) for L in (L1, L2, L3)]))
