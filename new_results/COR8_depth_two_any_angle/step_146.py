
import sys, pickle, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
R = D['R']
li4v = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve256.pkl','rb'))
bwords = list(itertools.product(range(4), repeat=4))
bidx = {w: i for i, w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(4) if l[j]]
def to256(S, e=0):
    v = [Fr(0)]*256
    for w, cf in S.items():
        coeff = cf.get(e)
        if coeff is None: continue
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for pp, c0 in terms:
            v[bidx[pp]] += coeff*c0
    return v
Rv = {E: to256(R, E) for E in (-1, 0, 1)}
li4v2, prodv, log2v = li4v['li4v'], li4v['prodv'], li4v['log2v']
dt = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/dt_syms_2kept.pkl','rb'))
import numpy as np
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
allv = list(li4v2.values()) + list(prodv.values()) + list(log2v.values()) + list(dt.values())
ra = rank_modp(allv)
print("rank all:", ra, " R_E in span:", all(rank_modp(allv + [Rv[E]]) == ra for E in (-1,0,1)))
