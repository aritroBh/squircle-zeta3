
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from itertools import product
from g2b_dbg import s05
from cor6_sym import sshuffle, rref

# weight-5 quotient T^5/D^5 over 2 letters
e0, e1 = (Fr(1), Fr(0)), (Fr(0), Fr(1))
B5 = list(product([e0, e1], repeat=5))
Drows = []
for a in (e0, e1):
    for bc in product([e0, e1], repeat=4):
        Drows.append([sshuffle({(a,): Fr(1)}, {bc: Fr(1)}).get(w, Fr(0)) for w in B5])
for ab in product([e0, e1], repeat=2):
    for cd in product([e0, e1], repeat=3):
        Drows.append([sshuffle({ab: Fr(1)}, {cd: Fr(1)}).get(w, Fr(0)) for w in B5])
Dmat, piv = rref([r[:] for r in Drows])
free = [j for j in range(32) if j not in piv]
print("dim D^5:", len(Dmat), " dim quotient:", len(free))

def qcoords(S):
    v = [S.get(w, Fr(0)) for w in B5]
    for i, pj in enumerate(piv):
        if v[pj]:
            q = v[pj]
            v = [a - q * b for a, b in zip(v, Dmat[i])]
    return tuple(v[j] for j in free)

def li5_sym_class(z):
    S = {(s05.vec_q(1 - z), s05.vec_q(z), s05.vec_q(z), s05.vec_q(z), s05.vec_q(z)): (Fr(-1),)}
    SE = s05.ft_expand_sym(S, 2)
    return qcoords({w: cf[0] for w, cf in SE.items()})

# target class: Omega_5(1/2)
S5 = s05.omega5_half()
tc = qcoords({w: cf[0] for w, cf in s05.ft_expand_sym(S5, 2).items()})
print("class(R4^(5)(1/2)):", [str(x) for x in tc])

args = [Fr(-2), Fr(1,3), Fr(2,3), Fr(1,2), Fr(1,9), Fr(-1,2), Fr(-1,3), Fr(1,4), Fr(3,4), Fr(-1,8), Fr(8,9), Fr(-1)]
cls = {}
for a in args:
    try:
        cls[a] = li5_sym_class(a)
    except Exception as ex:
        print("skip", a, ex)
print("Li5 classes:")
for a, v in cls.items():
    print(f"  {a}: {[str(x) for x in v]}")
# rank of Li5 span and membership
import sympy as sp
M = sp.Matrix([[sp.Rational(str(x)) for x in v] for v in cls.values()]).T
tv = sp.Matrix([sp.Rational(str(x)) for x in tc])
print("rank Li5 span:", M.rank())
aug = M.row_join(tv)
print("rank with target:", aug.rank(), "-> target in span:", M.rank() == aug.rank())
if M.rank() == aug.rank():
    sol = M.gauss_jordan_solve(tv)[0]
    for a, s in zip(cls.keys(), sol):
        if s: print(f"   {s} * Li5({a})")
