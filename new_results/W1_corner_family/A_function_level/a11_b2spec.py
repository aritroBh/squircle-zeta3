"""a11_b2spec.py -- 'formal specialisation' of the function-level Bloch group of an alphabet at c = 1/sqrt2.

B2 elements {f} for strictly admissible f (alphabet A); Lambda^2 image (1-f)^f (exact).  Function-level
relations = integer combinations with zero Lambda^2 image.  Point labels: the value f(1/sqrt2) modulo the
6-fold symmetry {x} = {1-1/x} = {1/(1-x)} = -{1/x} = -{1-x} = -{x/(x-1)} (labels carry the sign).
B2spec = (formal span of point labels) / (specialised function-level relations).
If rank B2spec = 1 then every point pair can be moved to a single generator using only function-level
relations + coincidences, i.e. a vanishing-pair realisation of ANY obstruction exists in this alphabet."""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
from sym2 import SymbolRing, b2
from admissible import enumerate_admissible
import flint
mp.mp.dps = 40
c = sp.Symbol('c')
R = SymbolRing([c])
roots = [Fr(0), Fr(1), Fr(-1), Fr(-1, 2), Fr(1, 2)]
BASE = [c, 1 - c, 1 + c, 1 + 2*c, 1 - 2*c]
for f in BASE: R.letter_index(f)
E = int(sys.argv[1]) if len(sys.argv) > 1 else 2
adm = []
for s, k, ex in enumerate_admissible(roots, E=E, K=3):
    f = sp.Integer(s)*sp.Integer(2)**k
    for r, e in zip(roots, ex):
        f *= (c - sp.Rational(r.numerator, r.denominator))**e
    adm.append(sp.cancel(f))
r2 = mp.sqrt(2)
def canon_label(v):
    """orbit representative + sign of {v}: choose the orbit element of smallest |.|+lexico, sign rel. to v"""
    orb = [(v, 1), (1 - 1/v, 1), (1/(1 - v), 1), (1/v, -1), (1 - v, -1), (v/(v - 1), -1)]
    best = min(orb, key=lambda t: (round(float(t[0]), 12)))
    return mp.nstr(best[0], 25), best[1]
labels = {}
rows = []
lam = []
for f in adm:
    v = sp.lambdify(c, f, 'mpmath')(1/r2)
    if abs(v) < mp.mpf(10)**-30 or abs(v - 1) < mp.mpf(10)**-30:
        continue
    lab, sg = canon_label(v)
    labels.setdefault(lab, len(labels))
    rows.append((str(f), lab, sg))
    lam.append(b2(R, f))
keys = sorted(set().union(*[set(d) for d in lam]))
nF = len(rows)
# matrix of Lambda^2 images: columns = f
A = flint.fmpq_mat(len(keys), nF)
for j, d in enumerate(lam):
    for i, k in enumerate(keys):
        v = d.get(k, 0)
        A[i, j] = flint.fmpq(v.numerator, v.denominator) if isinstance(v, Fr) else flint.fmpq(int(v))
# kernel of A (function-level relations)
rref, rank = A.rref()
piv = []
for i in range(rank):
    for j in range(nF):
        if rref[i, j] != 0:
            piv.append(j); break
free = [j for j in range(nF) if j not in piv]
print('admissible f (nondegenerate at point):', nF, ' Lambda^2-rank of B2:', rank, ' relations (kernel dim):', len(free))
# specialise relations to label space
L = len(labels)
specrows = []
for fj in free:
    vec = [flint.fmpq(0)]*nF
    vec[fj] = flint.fmpq(1)
    for i, pj in enumerate(piv):
        vec[pj] = -rref[i, fj]
    lv = [flint.fmpq(0)]*L
    for j in range(nF):
        if vec[j] != 0:
            lab, sg = rows[j][1], rows[j][2]
            lv[labels[lab]] += vec[j]*sg
    specrows.append(lv)
Ms = flint.fmpq_mat(len(specrows), L)
for i, rv in enumerate(specrows):
    for j in range(L): Ms[i, j] = rv[j]
rs = Ms.rref()[1]
print('point labels:', L, ' rank of specialised relations:', rs, '  => rank B2spec =', L - rs)
inv = {v: k for k, v in labels.items()}
print('labels:', [inv[i] for i in range(L)])
