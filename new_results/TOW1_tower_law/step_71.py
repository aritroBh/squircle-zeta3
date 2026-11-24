from towlib import *
# Build Gm_30^(k), Gm_21^(k) from the tower law (e-family) for k=3..14, compare affine law + slopes with WT3.
def Gm_from_towerlaw(key, k):
    """Gm as PF, using U_e tower values for depths 0..2d-4... need ALL depths for Gm!
       For degree 3 keys: depths 0,1,2 only (max 2d-4 = 2). So full tower known."""
    i, j = key; d = i+j
    maxdelta = 2*d-4
    ipos = [ZERO]
    for delta in range(0, maxdelta+1):
        r = 2*d-2-delta
        while len(ipos) < r+1: ipos.append(ZERO)
        ipos[r] = U_e(key, k, delta)
    hp = pf_zero(); hp[1] = ipos[1:] if ipos[0].is_zero() else ipos  # ensure index r-1 -> order r
    # rebuild properly: hp[1][r-1] = coeff of (v-i)^{-r}
    hp = pf_zero()
    hp[1] = [ipos[r] if r < len(ipos) else ZERO for r in range(1, 2*d-1)]
    hp[2] = [e.conj_i() for e in hp[1]]
    while hp[1] and hp[1][-1].is_zero(): hp[1].pop()
    while hp[2] and hp[2][-1].is_zero(): hp[2].pop()
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]

def Gm_actual(key, k):
    H = tower(k)[key]
    hp = pf_zero()
    maxr = max(H)
    hp[1] = [unq4(H.get(r,(0,0,0,0))) for r in range(1, maxr+1)]
    hp[2] = [e.conj_i() for e in hp[1]]
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]

# sanity: tower-law Gm == actual Gm for degree-3 keys
ok = True
for k in range(3, 15):
    for key in ((3,0),(2,1)):
        ok &= (pf_to_sympy(Gm_from_towerlaw(key,k)) == pf_to_sympy(Gm_actual(key,k)))
print('tower-law Gm == actual Gm, k=3..14:', ok)
# affine law + slopes
s30 = pf_from_sympy(sp.nsimplify(0)*0) if False else None
import sympy as sp
from fastred import v as fv
s30e = -(sp.Rational(1,12))*(1 - 2*sp.sqrt(2)*fv/(fv**2+1))
s30pf = pf_from_sympy(s30e)
s21pf = pf_scale(s30pf, Fr(-3))
ok2 = True
for k in range(4, 15):
    for key, s in ((3,0),s30pf), ((2,1),s21pf):
        lhs = pf_sub(Gm_from_towerlaw(key,k), Gm_from_towerlaw(key,3))
        rhs = pf_scale(s, Fr(k-3))
        if pf_to_sympy(lhs) != pf_to_sympy(rhs):
            ok2 = False; print('FAIL', key, k)
print('affine Gm law with WT3 slopes s30, s21: k=4..14:', ok2)
