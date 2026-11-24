import cor8r_sym as S
import sympy as sp
from fractions import Fraction
import itertools
cc = S.cc
A = 1/(1+cc); B = 2/(1+cc)
L = 4
values = [A, B]
unk = S._bad_words(L, values)
idx = {w:i for i,w in enumerate(unk)}
print("num unknowns:", len(unk))
alphabet = [sp.Integer(0)]+values
rows = []
for lu in range(1,L):
    lv = L-lu
    for u in itertools.product(alphabet, repeat=lu):
        su = S.symG(u)
        for v in itertools.product(alphabet, repeat=lv):
            sv = S.symG(v)
            prod = S.S_shuffle(su,sv)
            row = [Fraction(0)]*len(unk)
            rhs = dict(prod)
            ok = False
            for w in S._shuffle_words(u,v):
                if w in idx:
                    row[idx[w]] += 1; ok = True
                else:
                    sw = S.symG(w)
                    rhs = S.S_add(rhs, {ww:S.LP_neg(cf) for ww,cf in sw.items()})
            if ok:
                rows.append((u,v,row,rhs))
print("num rows:", len(rows))
# rank of coefficient matrix
M = sp.Matrix([r[2] for r in rows])
print("rank:", M.rank(), "unknowns:", len(unk))
# check consistency row by row symbolically: does the true assignment satisfy?
# true values: evaluate each unknown via a trustworthy numeric method? Instead, find rows that are linearly
# dependent but with different rhs. Compute RREF of augmented system over Q with rhs replaced by unique marker per row.
# Simpler: check whether each equation is a VALID shuffle identity for the true regularized values:
# evaluate G(w;1) numerically for the alphabet at c=0.31 via series? hard for bad words.
# Instead: find linear dependencies among coefficient rows and compare rhs.
import sympy as sp2
Mf = sp.Matrix([[Fraction(x) for x in r[2]] for r in rows])
ns = Mf.nullspace()
print("left nullity (row dependencies):", len(Mf.T.nullspace()))
