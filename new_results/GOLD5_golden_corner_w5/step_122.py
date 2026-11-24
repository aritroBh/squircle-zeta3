import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 30
P = G.golden_point('gol1')
F = P.F
ONE, ZERO = F.one(), F.zero()
def fval(x):
    return float(mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*mp.sqrt(5))
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
def li21_word(x, y):
    return (ZERO, F.inv(x), F.inv(F.mul(x, y)))
atoms = sorted({(x, y) for x in sols for y in sols
                if not F.eq(x, ONE) and any(F.eq(F.mul(x, y), s) for s in sols)
                and abs(fval(x)) <= 1 and abs(fval(F.mul(x, y))) < 1},
               key=lambda t: (float(abs(fval(t[0]))), float(abs(fval(t[1])))))
print(f"{len(atoms)} convergent Li21 atoms")
cl = {}
for x, y in atoms:
    S = G.make_S_J(P)(li21_word(x, y), ONE, ZERO)
    SE = G.ft_expand_sym_nf(S, 2)
    qc, _ = G.quotient_coords_nf(SE)
    cl[(str(x), str(y))] = qc
nz = {k: v for k, v in cl.items() if any(any(q) for q in v)}
print(f"{len(nz)} with nonzero Lie3 class")
for k, v in list(nz.items())[:15]:
    print(f"  Li21{k}: ({G.nf_str(v[0])}), ({G.nf_str(v[1])})")
# rank over Q (4-dim: 2 words x NF)
def q4(v):
    return (v[0][0], v[0][1], v[1][0], v[1][1])
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
from cor6_sym import rref
cols = []
for k, v in nz.items():
    cols.append(q4(v)); cols.append(q4((nfmul(SQ5v, v[0]), nfmul(SQ5v, v[1]))))
Mrref, piv = rref([list(r) for r in zip(*cols)])
print("Li21 Lie3-class span rank over Q:", len(piv))
json.dump({f"{k[0]}|{k[1]}": [[str(c) for c in q] for q in v] for k, v in cl.items()},
          open(GOLD + 'li21_classes.json', 'w'), indent=1)
