check = '''

import sympy as sp
from fastred import pf_to_sympy, pf_mul_DL, pf_fold, pf_scale, pf_add, DL1, DL2
v = sp.Symbol('v')

def Gm_of_key(i, j, k):
    H = tower(k)[(i,j)]
    hp = pf_zero()
    maxr = max(H) if H else 0
    ipos = [unq4(H.get(r,(0,0,0,0))) for r in range(1, maxr+1)]
    hp[1] = ipos
    hp[2] = [e.conj_i() for e in ipos]
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]

def diag_feed_reduction(p, k):
    """total feed from diagonal key (p,p) into child (p,p-1), as in the reducer."""
    Gm = Gm_of_key(p,p,k)
    direct = pf_scale(pf_mul_DL(Gm, DL2), Fr(-p))
    folded = pf_fold(pf_scale(pf_mul_DL(Gm, DL1), Fr(-p)))
    return pf_add(direct, folded)

def diag_feed_formula(p, k):
    """-p DL2 [Gm(v) - Gm(1/v)] as sympy."""
    Gm = Gm_of_key(p,p,k)
    G = pf_to_sympy(Gm)
    DL2e = 1/(v-AQ.sympy()) + 1/(v+BQ.sympy()) - 1/v
    return sp.expand(-p * DL2e * (G - G.subs(v, 1/v)))

for p in (2,3):
    for k in (2*p, 2*p+1):
        lhs = sp.expand(pf_to_sympy(diag_feed_reduction(p,k)))
        rhs = diag_feed_formula(p,k)
        print(f'p={p} k={k}: reduction-feed == -p DL2[Gm(v)-Gm(1/v)]:', sp.simplify(lhs-rhs)==0)
'''
open('towlib.py','a').write(check)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
