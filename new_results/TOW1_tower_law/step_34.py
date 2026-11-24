from towlib import *
import sympy as sp
from fastred import pf_to_sympy, pf_fold, pf_mul_DL, pf_scale, DL1, DL2, v as fv
def iszero(e):
    n, d = sp.fraction(sp.together(e))
    return sp.expand(n) == 0
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
aq, bq = AQ.sympy(), BQ.sympy()
DL2e = 1/(fv-aq) + 1/(fv+bq) - 1/fv
ok = True
for p in (2,3):
    for k in (2*p, 2*p+1, 2*p+2):
        Gm = Gm_of_key(p,p,k)
        direct = pf_scale(pf_mul_DL(Gm, DL2), Fr(-p))
        folded = pf_fold(pf_scale(pf_mul_DL(Gm, DL1), Fr(-p)))
        lhs = pf_to_sympy(pf_add(direct, folded))
        G = pf_to_sympy(Gm)
        rhs = -p*DL2e*(G - G.subs(fv,1/fv))
        r = iszero(lhs-rhs)
        ok &= r
        print(f'p={p} k={k}: child diag feed == -p DL2[Gm(v)-Gm(1/v)]: {r}', flush=True)
print('ALL:', ok)
