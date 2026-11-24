
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81_v2.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
import json
sa = json.load(open(C4+"sunit_args.json"))
# args with nonnegative c-exponent (no pole at c=0): check via sympy
import sympy as sp
c = sp.Symbol('c')
def c_exponent(s):
    ex = sp.sympify(s)
    num, den = sp.fraction(sp.together(ex))
    return sp.degree(sp.Poly(den, c)) - 0  # pole order at 0
polefree = []
for nm in li4v:
    ex = sp.sympify(nm)
    if ex.subs(c, 0).is_finite:
        polefree.append(nm)
print(len(polefree), "pole-free Li4 args:", polefree)
p = (1<<31)-1
def rank_modp(vecs):
    basis = []
    for v in vecs:
        v = np.array([x % p for x in v], dtype=np.int64)
        for b, pw in basis:
            if v[pw]: v = (v - int(v[pw])*b) % p
        nz = np.nonzero(v)[0]
        if len(nz):
            pw = int(nz[0]); basis.append(((v*pow(int(v[pw]), p-2, p)) % p, pw))
    return len(basis)
sub = [li4v[k] for k in polefree] + list(prodv.values())
ra = rank_modp(sub)
print("rank pole-free Li4 + products:", ra, " R_E in span:", all(rank_modp(sub + [Rv[E]]) == ra for E in (-1,0,1)))
