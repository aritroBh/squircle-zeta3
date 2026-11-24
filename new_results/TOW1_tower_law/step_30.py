# neutralize the auto-running check block in towlib.py
src = open('towlib.py').read()
i = src.find("\nimport sympy as sp\nfrom fastred import pf_to_sympy")
if i > 0:
    src = src[:i] + "\n"   # drop the check block; keep library part
    open('towlib.py','w').write(src)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
import sympy as sp
from fastred import pf_to_sympy, pf_mul_DL, pf_fold, pf_scale, DL1, DL2
v = sp.Symbol('v')
aq, bq = AQ.sympy(), BQ.sympy()
DL1e = 1/(v-bq) + 1/(v+aq) - 1/v
DL2e = 1/(v-aq) + 1/(v+bq) - 1/v
print('DL1(1/v)/v^2 == -DL2(v):', sp.simplify(sp.expand(DL1e.subs(v,1/v)/v**2 + DL2e))==0)

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

p,k = 2,4
Gm = Gm_of_key(p,p,k)
R = pf_scale(pf_mul_DL(Gm, DL1), Fr(-p))
Re = pf_to_sympy(R)
print('pf_fold(R) == R(1/v)/v^2:', sp.simplify(sp.expand(pf_to_sympy(pf_fold(R)) - Re.subs(v,1/v)/v**2))==0)
direct = pf_scale(pf_mul_DL(Gm, DL2), Fr(-p))
de = pf_to_sympy(direct)
GmE = pf_to_sympy(Gm)
formula = -p*DL2e*(GmE - GmE.subs(v,1/v))
lhs = sp.expand(de + Re.subs(v,1/v)/v**2)
print('reduction diag feed == -p DL2[Gm(v)-Gm(1/v)]:', sp.simplify(sp.expand(lhs - formula))==0)
