import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 40
r5 = mp.sqrt(5)

P = G.golden_point('gol1')
F = P.F
ONE = F.one()
def fval(x):
    return mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5
def vec_or_none(el):
    try: return P.vec(el)
    except ValueError: return None
sols = [ONE, F.neg(ONE)]
for tup in itertools.product(range(-9, 10), repeat=2):
    if tup == (0,0): continue
    v = ONE
    for e_, g in zip(tup, G.BASIS):
        if e_: v = F.mul(v, F.pow(g, e_))
    for sgn in (1, -1):
        x = v if sgn == 1 else F.neg(v)
        if F.eq(x, ONE): continue
        om = F.sub(ONE, x)
        if all(t == 0 for t in om): continue
        if vec_or_none(om) is not None and not any(F.eq(x, y) for y in sols):
            sols.append(x)
zs = [s for s in sols if not F.eq(s, ONE)]
print(f"{len(zs)} Li5 arguments")
li5cl = {}
for z in zs:
    cl = G.quotient5(G.ft_expand_sym_nf(G.li5_sym_nf(P, z), 2))
    li5cl[z] = cl
nz = {z: cl for z, cl in li5cl.items() if any(any(q) for q in cl)}
print(f"{len(nz)} with nonzero class")
import json
json.dump({str(z): [[str(c) for c in q] for q in cl] for z, cl in li5cl.items()},
          open(GOLD + 'li5_classes_golden.json', 'w'), indent=1)
# rank of the span over Q (12-dim: per z cols v_z, sqrt5*v_z)
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
cols = []
names = []
for z, cl in li5cl.items():
    cols.append(q12(cl)); names.append((z, 'a'))
    cols.append(q12(tuple(nfmul(SQ5v, q) for q in cl))); names.append((z, 'b'))
M = [list(c) for c in cols]
R = G.s01 and None
from cor6_sym import rref
Mrref, piv = rref([list(r) for r in zip(*M)])  # rows = Q-coords
# rank:
print("Li5-class span rank over Q:", len(piv))
json.dump({'zs': [str(z) for z in zs]}, open(GOLD + 'li5_args.json', 'w'))
