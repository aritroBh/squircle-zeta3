
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from itertools import product
from g2b_dbg import s05
from cor6_sym import sshuffle, rref
import sympy as sp

B5 = list(product(range(2), repeat=5))
Drows = []
for a in range(2):
    for bc in product(range(2), repeat=4):
        Drows.append([sshuffle({(a,): Fr(1)}, {bc: Fr(1)}).get(w, Fr(0)) for w in B5])
for ab in product(range(2), repeat=2):
    for cd in product(range(2), repeat=3):
        Drows.append([sshuffle({ab: Fr(1)}, {cd: Fr(1)}).get(w, Fr(0)) for w in B5])
Dmat, piv = rref([r[:] for r in Drows])
free = [j for j in range(32) if j not in piv]

def qcoords(S):
    v = [S.get(w, Fr(0)) for w in B5]
    for i, pj in enumerate(piv):
        if v[pj]:
            q = v[pj]
            v = [a - q * b for a, b in zip(v, Dmat[i])]
    return tuple(v[j] for j in free)

def cls_of_expanded(SE):
    return qcoords({w: cf[0] for w, cf in SE.items()})

def li5_sym_class(z):
    S = {(s05.vec_q(1 - z), s05.vec_q(z), s05.vec_q(z), s05.vec_q(z), s05.vec_q(z)): (Fr(-1),)}
    return cls_of_expanded(s05.ft_expand_sym(S, 2))

# Li32 atom classes
def li32_class(x, y):
    S = s05.li32_sym(x, y)
    return cls_of_expanded(s05.ft_expand_sym({w: (cf,) for w, cf in S.items()}, 2))

S5 = s05.omega5_half()
tc = qcoords({w: cf[0] for w, cf in s05.ft_expand_sym(S5, 2).items()})
print("class(R4^(5)):", [str(x) for x in tc])

atoms_N5 = [(Fr(1,4), Fr(4,3), Fr(5,4)), (Fr(8,9), Fr(1), Fr(1,3))]
n5cls = [Fr(0)]*6
for x, y, d in atoms_N5:
    cl = li32_class(x, y)
    print(f"class Li32({x},{y}): {[str(x_) for x_ in cl]}")
    n5cls = [a + d*Fr(b) for a, b in zip(n5cls, cl)]
print("class(N5):", [str(x) for x in n5cls])
rest = [a - b for a, b in zip(tc, n5cls)]
print("class(R4^(5) - N5):", [str(x) for x in rest])

args = [Fr(-2), Fr(1,3), Fr(2,3), Fr(1,2), Fr(1,9), Fr(-1,2), Fr(-1,3), Fr(1,4), Fr(3,4), Fr(-1,8), Fr(8,9), Fr(-1)]
cls = {a: li5_sym_class(a) for a in args}
M = sp.Matrix([[sp.Rational(str(x)) for x in v] for v in cls.values()]).T
rv = sp.Matrix([sp.Rational(str(x)) for x in rest])
print("Li5 span rank:", M.rank(), " rank with (R4-N5):", M.row_join(rv).rank())
if M.rank() == M.row_join(rv).rank():
    sol = M.gauss_jordan_solve(rv)[0]
    for a, s in zip(cls.keys(), sol):
        if s: print(f"   {s} * Li5({a})")
