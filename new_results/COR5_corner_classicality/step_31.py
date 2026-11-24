import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *
from itertools import product as iproduct

e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
letters = [e1, e2]
def words_of_len(n, lets): return list(iproduct(lets, repeat=n))
W1 = words_of_len(1, letters); W2 = words_of_len(2, letters)
W3 = words_of_len(3, letters); W4 = words_of_len(4, letters)

# decomposable subspace D^4 = span of all shuffle products u ш v, |u|+|v|=4
D_rows = []
for u in W1:
    for v in W3:
        D_rows.append(sshuffle({u: Fr(1)}, {v: Fr(1)}))
for u in W2:
    for v in W2:
        D_rows.append(sshuffle({u: Fr(1)}, {v: Fr(1)}))
for u in W3:
    for u_ in W1:
        pass  # (3,1) same as (1,3) by commutativity of shuffle
def vec_of(S): return tuple(S.get(w, Fr(0)) for w in W4)

def rref(rows):
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
    return rows[:R], piv

Dmat, pivD = rref([vec_of(s) for s in D_rows])
print("dim D^4 (decomposables):", len(Dmat), " (expect 13); quotient dim:", 16 - len(Dmat))
free_cols = [j for j in range(16) if j not in pivD]
print("quotient basis words:", [W4[j] for j in free_cols])

# quotient coordinates of a symbol S: reduce vec against Dmat, keep free columns
def quotient_coords(S):
    v = list(vec_of(S))
    for i, pj in enumerate(pivD):
        if v[pj]:
            q = v[pj]
            v = [a - q*b for a, b in zip(v, Dmat[i])]
    return tuple(v[j] for j in free_cols)

# test: Li4 symbol nonzero in quotient
rv = e1; om = e2
S_li4 = li_sym(4, rv, om)
print("S(Li4(x)) mod products:", quotient_coords(S_li4))
# shuffle product zero in quotient
S_sh = sshuffle({W1[0]: Fr(1)}, {W3[3]: Fr(1)})
print("shuffle mod products:", quotient_coords(S_sh))
