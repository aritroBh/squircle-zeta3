# route 2: assemble D~_{n,m} = sum_j 2 e_{m,n-1-j} I_{n-1-j,j}/j! from the (exactly certified) I-reductions; non-zeta symbols must cancel formally
from common import *
import sympy as sp
d89 = pickle.load(open(T+'I_reductions_w89.pkl', 'rb')); d9 = pickle.load(open(T+'I_reductions_w9.pkl', 'rb'))
basis = {8: d89['basis8'], 9: d89['basis9']}; reds = {8: d89['reds'], 9: d9}
Dtarget = pickle.load(open(T+'Dtilde_relations.pkl', 'rb'))
def is_zeta(nm): return all(f[0] == 'z' and not f.startswith('zbar') for f in nm.split('*'))
z2k = lambda k: sp.zeta(2*k) / sp.zeta(2)**k      # zeta(2k)/zeta(2)^k  (rational)
def zeta_monomial_to_target(nm, W):
    """express a z-monomial of the basis in the 'Dtarget' basis ZB[W]: only change is z2^k -> rational * zeta(2k) (W even, k=W/2)"""
    fs = nm.split('*'); ns = tuple(sorted(int(f[1:]) for f in fs))
    if all(n == 2 for n in ns):
        return (W,), F(str(sp.nsimplify(sp.zeta(W)/sp.zeta(2)**(W//2)))) ** -1  # z2^k = (z2^k/zk) zk
    return ns, F(1)
for W in (8, 9):
    for m in range(W//2):
        acc = {}
        for j in range(1, W-2*m):
            r = W-1-j; e = e_table(m, r)
            if e == 0: continue
            c = 2*e/fact(j)
            for cf, nm in zip(reds[W][(r, j)], basis[W]):
                if cf: acc[nm] = acc.get(nm, F(0)) + c*F(cf)
        nonz = {nm: v for nm, v in acc.items() if v and not is_zeta(nm)}
        zet = {}
        for nm, v in acc.items():
            if v and is_zeta(nm):
                ns, fac = zeta_monomial_to_target(nm, W); zet[ns] = zet.get(ns, F(0)) + v*fac
        tgt = {tuple(sorted(ns)): F(c) for ns, c in zip(ZB[W], Dtarget[f"{W},{m}"]) if F(c)}
        zet = {k: v for k, v in zet.items() if v}
        print(f"D~_{{{W},{m}}}: non-zeta symbols left: {len(nonz)}; zeta part == stored target: {zet == tgt}")
