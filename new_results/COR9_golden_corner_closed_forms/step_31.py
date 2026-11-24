import time, itertools
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_lib import QuadField, QuadPoint, li4_coords, rankQ

QF5 = QuadField(5)
P = QuadPoint(QF5, "1/phi", (Fr(-1,2), Fr(1,2)), primes=[2])
F = P.F
print("gens:", [(kind, l) for _, kind, l in P.gensE])
gens = [g[0] for g in P.gensE]
for el, kind, l in P.gensE:
    print("  gen", kind, el, "N =", F.norm(el))

def sunit_args(P, K=9):
    F = P.F
    gens = [g[0] for g in P.gensE]
    args = []
    for tup in itertools.product(range(-K, K+1), repeat=len(gens)):
        v = F.one()
        for e_, g in zip(tup, gens):
            if e_: v = F.mul(v, F.pow(g, e_))
        for sgn in (1, -1):
            x = v if sgn == 1 else F.neg(v)
            if F.eq(x, F.one()): continue
            om = F.sub(F.one(), x)
            if all(t == 0 for t in om): continue
            if P.vec(om) is not None and not any(F.eq(x, y) for y in args):
                args.append(x)
    return args

t0 = time.time()
args = sunit_args(P, K=9)
print(f"\n{len(args)} {2}-unit Li4 arguments in Q(sqrt5) ({time.time()-t0:.0f}s):")
import mpmath as mp
mp.mp.dps = 30
r5 = mp.sqrt(5)
def fval(x): return float(mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5)
li4cl = []
for x in sorted(args, key=fval):
    v = li4_coords(P, x)
    li4cl.append((x, v))
    print(f"  x = {str(x):28s} ~ {fval(x):+12.6f}   class {v}")
print("rank of Li4 classes:", rankQ([v for _, v in li4cl]))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_31: Q(sqrt5) S-unit Li4 args enumerated: {len(args)}, class rank {rankQ([v for _, v in li4cl])}.\n")
