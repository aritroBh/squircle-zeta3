import sys, itertools
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
from gold5_lib import golden_point, PHI, TWO
import mpmath as mp
mp.mp.dps = 30
r5 = mp.sqrt(5)
P = golden_point('gol1')
F = P.F
def fval(x):
    return mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5
gens = [TWO, PHI]
args = []
K = 9
for tup in itertools.product(range(-K, K+1), repeat=2):
    if tup == (0,0): continue
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
print(f"{len(args)} exceptional {{2}}-units in Q(sqrt5):")
sols = []
for x in sorted(args, key=fval):
    vx = P.vec(x); vo = P.vec(F.sub(F.one(), x))
    sols.append((x, vx, vo))
    print(f"  x={str(x):24s} ~ {float(fval(x)):+12.8f}  v(x)={tuple(str(t) for t in vx)} v(1-x)={tuple(str(t) for t in vo)}")
