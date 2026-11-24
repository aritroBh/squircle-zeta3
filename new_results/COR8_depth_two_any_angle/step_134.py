
import pickle
from fractions import Fraction as Fr
import sympy as sp
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve256.pkl','rb'))
li4v, prodv, log2v, Rv = D['li4v'], D['prodv'], D['log2v'], D['Rv']
allv = list(li4v.values()) + list(prodv.values()) + list(log2v.values())
M = sp.Matrix(256, len(allv), lambda i, j: sp.Rational(allv[j][i]))
# solve M x = R_E per exponent
import itertools
bwords = list(itertools.product(range(4), repeat=4))
for E in (-1, 0, 1):
    b = sp.Matrix(256, 1, lambda i, j: sp.Rational(Rv[E][i]))
    r = M.row_join(b).rref()[0]
    inc = [i for i in range(256) if all(r[i, j] == 0 for j in range(len(allv))) and r[i, len(allv)] != 0]
    print(f"E={E}: {len(inc)} unsatisfied components")
    for i in inc[:20]:
        bw = bwords[i]
        print("   word", ''.join('2cmp'[x] for x in bw), "target", sp.nsimplify(Rv[E][i]))
