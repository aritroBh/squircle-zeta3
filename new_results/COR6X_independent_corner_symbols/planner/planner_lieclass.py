# Planner check: are the paper's two Li4 combinations at each classical {2,3}-corner in the same Lie class (symbol modulo shuffles)?
# Letters modulo torsion: [2^a 3^b] = (a, b). Symbol of Li4(x) = -[1-x] (x) [x] (x) [x] (x) [x].
import itertools
from fractions import Fraction as Fr
from sympy import Matrix, factorint, Rational
def vec(q):
    q = Fr(q); f = {}
    for p, k in factorint(abs(q.numerator)).items(): f[p] = f.get(p, 0) + k
    for p, k in factorint(q.denominator).items(): f[p] = f.get(p, 0) - k
    assert set(f) <= {2, 3}, q
    return (f.get(2, 0), f.get(3, 0))
def sym_li4(x):
    a, b = vec(1 - Fr(x)), vec(x)
    T = {}
    for w in itertools.product(range(2), repeat=4):
        T[w] = -a[w[0]] * b[w[1]] * b[w[2]] * b[w[3]]
    return T
def comb(terms):
    T = {w: Fr(0) for w in itertools.product(range(2), repeat=4)}
    for c, x in terms:
        for w, v in sym_li4(x).items(): T[w] += Fr(c) * v
    return T
def shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + s for s in shuffle(u[1:], v)] + [(v[0],) + s for s in shuffle(u, v[1:])]
words = list(itertools.product(range(2), repeat=4)); idx = {w: i for i, w in enumerate(words)}
gens = []
for k in (1, 2):
    for u in itertools.product(range(2), repeat=k):
        for v in itertools.product(range(2), repeat=4 - k):
            row = [0] * 16
            for s in shuffle(u, v): row[idx[s]] += 1
            gens.append(row)
Sh = Matrix(gens); r0 = Sh.rank(); print("shuffle subspace rank:", r0, "(quotient dim", 16 - r0, ")")
def same_class(A, B):
    d = [Rational(A[w] - B[w]) for w in words]
    return Matrix(gens + [d]).rank() == r0
half_paper = [(-21, -2), (Fr(13, 2), Fr(1, 3)), (Fr(13, 2), Fr(2, 3))]
half_numer = [(-21, -2), (Fr(-65, 8), Fr(1, 3)), (Fr(13, 2), Fr(2, 3)), (Fr(137, 6), Fr(1, 2)), (Fr(39, 64), Fr(1, 9))]
third_paper = [(Fr(-11, 81), -8), (-2, -3), (Fr(20, 3), Fr(2, 3))]
third_numer = [(Fr(-22, 3), -2), (Fr(-17, 2), Fr(1, 3)), (Fr(20, 3), Fr(2, 3)), (Fr(110, 9), Fr(1, 2)), (Fr(25, 48), Fr(1, 9))]
print("Omega(1/2): paper's 3-term Lie class == Li4 part of the 1000-digit closed form:", same_class(comb(half_paper), comb(half_numer)))
print("Omega(1/3): paper's 3-term Lie class == Li4 part of the 1000-digit W4 closed form:", same_class(comb(third_paper), comb(third_numer)))
print("sanity: Li4(1/2) alone is a nonzero class:", not same_class(comb([(1, Fr(1, 2))]), comb([])))
