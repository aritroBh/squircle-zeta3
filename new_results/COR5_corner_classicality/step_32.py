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

def expand_word(w):
    """word over V (letter vectors in Q^2) -> dict: word over basis letters -> Fr (multilinear)."""
    terms = [((), Fr(1))]
    for l in w:
        terms = [(p + (e1,), c * l[0]) for p, c in terms if l[0]] + \
                [(p + (e2,), c * l[1]) for p, c in terms if l[1]]
    out = {}
    for p, c in terms:
        out[p] = out.get(p, 0) + c
    return {p: c for p, c in out.items() if c}

def expand_sym(S):
    out = {}
    for w, cf in S.items():
        for p, c in expand_word(w).items():
            out[p] = out.get(p, 0) + cf * c
    return {p: c for p, c in out.items() if c}

# rebuild decomposables over basis words
D_rows = []
for u in W1:
    for v in W3: D_rows.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)})))
for u in W2:
    for v in W2: D_rows.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)})))
def vec_of(S): return tuple(S.get(w, Fr(0)) for w in W4)
def rref(rows):
    rows = [list(r) for r in rows]; R = 0; piv = []
    for j in range(len(rows[0])):
        p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
        if p is None: continue
        rows[R], rows[p] = rows[p], rows[R]
        q = rows[R][j]; rows[R] = [x/q for x in rows[R]]
        for i in range(len(rows)):
            if i != R and rows[i][j]:
                rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
        piv.append(j); R += 1
        if R == len(rows): break
    return rows[:R], piv
Dmat, pivD = rref([vec_of(s) for s in D_rows])
free_cols = [j for j in range(16) if j not in pivD]
print("dim D^4:", len(Dmat), "; quotient dim:", 16 - len(Dmat))
def quotient_coords(S):
    v = list(vec_of(expand_sym(S)))
    for i, pj in enumerate(pivD):
        if v[pj]:
            q = v[pj]
            v = [a - q*b for a, b in zip(v, Dmat[i])]
    assert all(v[j] == 0 for j in range(16) if j not in free_cols), "residual outside quotient"
    return tuple(v[j] for j in free_cols)

S_li4 = li_sym(4, e1, e2)
print("S(Li4(x)) mod products:", quotient_coords(S_li4), " (basis words:", [ "".join('1' if l==e1 else '2' for l in W4[j]) for j in free_cols], ")")
S_sh = sshuffle({W1[0]: Fr(1)}, {W3[3]: Fr(1)})
print("shuffle mod products:", quotient_coords(S_sh))
