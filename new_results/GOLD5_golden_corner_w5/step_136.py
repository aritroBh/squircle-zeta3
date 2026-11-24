import sys, itertools, json, pickle
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
from cor6_sym import rref
import mpmath as mp
mp.mp.dps = 30
D = pickle.load(open(GOLD + 'w3_content_gol1.pkl', 'rb'))
w3_om, w3_atoms, w3_li5, W3M = D['w3_om'], D['w3_atoms'], D['w3_li5'], D['W3M']
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v)
          for k, v in json.load(open(GOLD + 'atom_classes_golden.json')).items()}
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
keys = list(w3_atoms.keys())
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
nzsmall = [z for z in li5cl if abs(fval(z)) <= 1]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
cl_om1 = ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)),
          (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21)))
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
def pack_atom(k):
    nfv = [(V[k][0], Fr(0)), (V[k][1], Fr(0))]
    nfv += [pad(q) for q in atomcl[k]]
    nfv += [pad(q) for q in w3_atoms[k]]
    return nfv
def pack_li5(z):
    nfv = [(Fr(0), Fr(0)), (Fr(0), Fr(0))]
    nfv += [pad(q) for q in li5cl[z]]
    nfv += [pad(q) for q in w3_li5[z]]
    return nfv
Tom = [(Fr(-72), Fr(0)), (Fr(72), Fr(-54))] + [pad(q) for q in cl_om1] + [pad(q) for q in w3_om]
def qvec(nfv):
    out = []
    for a in nfv: out += [a[0], a[1]]
    return tuple(out)
Tq = qvec(Tom)
cols = []; names = []
for k in keys:
    nfv = pack_atom(k)
    cols.append(qvec(nfv)); names.append(('atom', k, 0))
    cols.append(qvec(tuple(nfmul(SQ5v, a) for a in nfv))); names.append(('atom', k, 1))
for z in nzsmall:
    nfv = pack_li5(z)
    cols.append(qvec(nfv)); names.append(('li5', z, 0))
    cols.append(qvec(tuple(nfmul(SQ5v, a) for a in nfv))); names.append(('li5', z, 1))
M = [list(c) for c in cols]
R1 = rref([list(r) for r in zip(*M)])[1]
R2 = rref([list(r) for r in zip(*(M + [list(Tq)]))])[1]
print(f"gol1: constraint dim {len(Tq)}, #cols {len(cols)}, rank {len(R1)}, rank+target {len(R2)}, solvable: {len(R1)==len(R2)}")
