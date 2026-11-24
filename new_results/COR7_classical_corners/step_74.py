import sys, itertools
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
from cor6_sym import W4Quotient, sadd, lneg

def sunit_li4_args(P, K=10):
    """x = +- prod gens^e (e in [-K,K]) with x != 0,1 and 1-x also in the group. Exact via P.vec."""
    F = P.F
    gens = [g[0] for g in P.gensE]
    n = len(gens)
    args = []
    for tup in itertools.product(range(-K, K+1), repeat=n):
        v = F.one()
        for e_, g in zip(tup, gens):
            if e_:
                v = F.mul(v, F.pow(g, e_))
        # x = v or -v; need x != 1 and vec(1-x)
        for sgn in (1, -1):
            x = v if sgn == 1 else F.neg(v)
            if F.eq(x, F.one()):
                continue
            om = F.sub(F.one(), x)
            if all(t == 0 for t in om):
                continue
            if P.vec(om) is not None:
                args.append((x, tup, sgn))
    # dedupe exact duplicates
    uniq = []
    for x, tup, sgn in args:
        if not any(F.eq(x, y) for y, _, _ in uniq):
            uniq.append((x, tup, sgn))
    return uniq

t_out = {}
for D, c0s in [(5, [(Fr(-1,2), Fr(1,2)), (Fr(-2), Fr(1))]), (2, [(Fr(-1), Fr(1)), (Fr(3), Fr(-2)), (Fr(0), Fr(1,2))])]:
    QF = QuadField(D)
    P = QuadPoint(QF, f"D{D}", c0s[0], primes=[2])
    args = sunit_li4_args(P, K=9)
    print(f"Q(sqrt{D}): {len(args)} S-unit Li4 arguments (x and 1-x in <gens>):")
    for x, tup, sgn in args:
        print("   x =", x, "=", ("-" if sgn<0 else "")+"*".join(f"g{i}^{e}" for i,e in enumerate(tup) if e), flush=True)
    t_out[D] = args
