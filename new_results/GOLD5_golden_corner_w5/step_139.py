import sys, itertools, json, pickle
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
from cor6_sym import rref
import mpmath as mp
mp.mp.dps = 30
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v)
          for k, v in json.load(open(GOLD + 'atom_classes_golden.json')).items()}
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
nzsmall = [z for z in li5cl if abs(fval(z)) <= 1]
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
def qvec(nfv):
    out = []
    for a in nfv: out += [a[0], a[1]]
    return tuple(out)
for tag in ('gol1', 'gol3'):
    D = pickle.load(open(GOLD + f'w3_content_{tag}.pkl', 'rb'))
    w3_li5 = D['w3_li5']
    # li5 column rank on class+w3 (12+16 = 28 Q-dims)
    cols = []
    for z in nzsmall:
        nfv = [pad(q) for q in li5cl[z]] + [pad(q) for q in w3_li5[z]]
        cols.append(qvec(nfv))
        cols.append(qvec(tuple(nfmul(SQ5v, a) for a in nfv)))
    M = [list(c) for c in cols]
    R = rref([list(r) for r in zip(*M)])[1]
    print(f"{tag}: li5-column rank on class+w3 (28 dims): {len(R)} ({len(nzsmall)} args)")
