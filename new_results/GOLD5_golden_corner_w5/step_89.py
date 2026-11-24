import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 30
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
small = [z for z, cl in li5cl.items() if abs(fval(z)) <= 1]
nzsmall = [z for z in small if any(any(q) for q in li5cl[z])]
print(f"{len(small)} args with |x|<=1: {[str(z) for z in small]}")
print(f"  of which nonzero class: {len(nzsmall)}")
for z in nzsmall:
    print(f"    z={str(z):24s} ~ {fval(z):+.6f}: {[G.nf_str(q) for q in li5cl[z]]}")
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]
cols = []
for z in nzsmall: cols += cols_z(z)
from cor6_sym import rref
Mrref, piv = rref([list(r) for r in zip(*cols)])
print("small-arg class span rank over Q:", len(piv))
