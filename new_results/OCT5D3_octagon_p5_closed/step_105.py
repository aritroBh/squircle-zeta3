import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from fractions import Fraction as F
import sympy as sp

names = pickle.load(open("basis_keep.pkl","rb"))['names']
rel1500 = {k: tuple(F(x) for x in v) for k, v in pickle.load(open("rel1500_frac.pkl","rb")).items()}
rel1100 = {k: tuple(F(x) for x in v) for k, v in pickle.load(open("rel1100_frac.pkl","rb")).items()}
def vec102(rel):
    v = [F(0)]*102
    for i, nm in enumerate(names):
        if nm in rel:
            v[2*i], v[2*i+1] = rel[nm]
    return v
d = [vec102(rel1500)[j] - vec102(rel1100)[j] for j in range(102)]
ak = pickle.load(open("atom_kernel.pkl","rb"))
rels = [v for h, v in ak['rels'] if h < 2e6]
print("true kernel relations used:", len(rels))
M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in rels])
M2 = M.row_insert(len(rels), sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in d]]))
print("rank(kernel) =", M.rank(), " rank(kernel + diff) =", M2.rank())
print("rel1500 - rel1100 lies in the verified atom-kernel:", M.rank() == M2.rank())
