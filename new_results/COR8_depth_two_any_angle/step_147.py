
import pickle
from fractions import Fraction as Fr
import sympy as sp
import itertools
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve256.pkl','rb'))
li4v, prodv, log2v = D['li4v'], D['prodv'], D['log2v']
dt = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/dt_syms_2kept.pkl','rb'))
Rv = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/Rv_2kept.pkl','rb'))
allv = list(li4v.values()) + list(prodv.values()) + list(log2v.values()) + list(dt.values())
M = sp.Matrix(256, len(allv), lambda i, j: sp.Rational(allv[j][i]))
ln = M.T.nullspace()
print("left nullity:", len(ln))
bwords = list(itertools.product(range(4), repeat=4))
for E in (-1,):
    bad = []
    for v in ln:
        val = sum(v[i]*Rv[E][i] for i in range(256))
        if val != 0:
            bad.append((v, sp.nsimplify(val)))
    print(f"E={E}: {len(bad)} violated")
    for v, val in bad[:10]:
        sup = [(''.join('2cmp'[x] for x in bwords[i]), sp.nsimplify(v[i])) for i in range(256) if v[i]]
        print("   val", val, ":", sup[:14])
