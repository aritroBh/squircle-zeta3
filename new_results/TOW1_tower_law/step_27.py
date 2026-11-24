tests = '''

def Gm_of_key(i, j, k):
    """recompute Gm_ij^(k) as PF from the tower data + eval_b: Gm = (1-eval_b) Int HigherPart(W)."""
    H = tower(k)[(i,j)]
    from fastred import pf_zero, pf_integrate_higher, pf_eval_at
    hp = pf_zero()
    maxr = max(H) if H else 0
    ipos = [unq4(H.get(r,(0,0,0,0))) for r in range(1, maxr+1)]
    hp[1] = ipos
    hp[2] = [e.conj_i() for e in ipos]
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]

def pf_subP(P, Q):
    from fastred import pf_add, pf_neg
    return pf_add(P, pf_neg(Q))

def diag_feed_rhs(p, k):
    """-p * DL2 * [Gm_pp(v) - Gm_pp(1/v)]  as PF."""
    from fastred import pf_mul_DL, pf_fold, pf_scale, pf_add
    Gm = Gm_of_key(p,p,k)
    F = pf_fold(Gm)   # Gm(1/v)/v^2 -- NOT what we want; we want Gm(1/v) plain.
    return None
'''
# pf_fold includes the /v^2; I need plain v->1/v substitution. Let me write pf_subst_inv.
subst = '''

def pf_inv_subst(P):
    """R(v) -> R(1/v) (NO /v^2), for PF with zero poly part; may create pole at 0 / poly part? 
       For our use Gm vanishes at infinity (O(1/v)), so R(1/v) has zero poly part."""
    real, ipos, ineg, const = P
    # use pf_fold-like logic but without /v^2: R(1/v) = F(R)*v^2 ... instead build directly:
    # simple pole c/(v-r) -> c/(1/v - r) = -c/r + (c/r^2)/(1/r - v)*... handle via sympy per-term? 
    # Easier: R(1/v) = F(R)(v) * v^2 where F = pf_fold. v^2 * PF: multiply.
    from fastred import pf_fold
    F = pf_fold(P)   # = R(1/v)/v^2
    return pf_mul_v2(F)
'''
open('/tmp/t.py','w').write('')
print('need v^2 multiplication helper; doing it with sympy conversion instead for checking')
