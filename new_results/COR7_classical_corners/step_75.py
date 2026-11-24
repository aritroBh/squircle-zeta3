import sys, itertools
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
from cor6_pipe import W4TERMS
from cor6_sym import W4Quotient, lneg
from cor6_ft import ft_coords

Q2 = W4Quotient(2)

def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})

def sunit_li4_args(P, K=9):
    F = P.F
    gens = [g[0] for g in P.gensE]
    args = []
    for tup in itertools.product(range(-K, K+1), repeat=len(gens)):
        v = F.one()
        for e_, g in zip(tup, gens):
            if e_:
                v = F.mul(v, F.pow(g, e_))
        for sgn in (1, -1):
            x = v if sgn == 1 else F.neg(v)
            if F.eq(x, F.one()):
                continue
            om = F.sub(F.one(), x)
            if P.vec(om) is not None:
                if not any(F.eq(x, y) for y in args):
                    args.append(x)
    return args

def rankQ(vecs):
    M = [list(v) for v in vecs]
    R = 0
    for c in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][c]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][c]:
                M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
        R += 1
    return R

CORNERS = [
    (5, "1/phi (golden)", (Fr(-1,2), Fr(1,2))),
    (5, "1/phi^3", (Fr(-2), Fr(1))),
    (2, "sqrt2-1", (Fr(-1), Fr(1))),
    (2, "(sqrt2-1)^2", (Fr(3), Fr(-2))),
    (2, "1/sqrt2 (square)", (Fr(0), Fr(1,2))),
]
store = {}
for D, nm, c0 in CORNERS:
    QF = QuadField(D)
    P = QuadPoint(QF, nm, c0, primes=[2])
    args = sunit_li4_args(P, K=8)
    li4 = [(x, li4_coords(P, x)) for x in args]
    rLi4 = rankQ([v for _, v in li4])
    SOm = P.omega_symbol(W4TERMS)
    cls = ft_coords(Q2, SOm, 2)
    # membership of each component in Q-span of li4
    span = [v for _, v in li4]
    in_span = []
    for comp in cls:
        in_span.append(rankQ(span + [list(comp)]) == rankQ(span))
    print(f"Q(sqrt{D}) {nm}: rank Li4 classes = {rLi4} (of {len(args)} args); class(Omega) components:")
    for j, comp in enumerate(cls):
        print(f"    theta^{j}: {tuple(str(x) for x in comp)}   in Q-span Li4: {in_span[j]}")
    store[(D, nm)] = (cls, li4)
import pickle
pickle.dump({k: (v[0], [(str(x_ ), tuple(str(c) for c in v_)) for x_, v_ in v[1]]) for k, v in store.items()}, open('lie_classes.pkl','wb'))
