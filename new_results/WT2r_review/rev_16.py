
import sys, time
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fastred import v
from fractions import Fraction as Fr
import buildintegrand
S2=sp.sqrt(2)
Rw2 = 1/(2*v) + 1/(2*(v-(1+S2))) - 1/(v+(1+S2)) - 1/(2*(v+(S2-1)))
# Rw2 denominators: v, v-a, v+a, v+b -- all real; no (v^2+1) factor
den = sp.fraction(sp.together(Rw2))[1]
print("8. Rw2 denominator factors:", sp.factor(den), "-> no v^2+1: +-i poles NONE;",
      " and (1-eps)W_10 = eps^2 Rw2 (rev_15 item 7) => W_10(eps) has no +-i poles for all eps")

# k=13 stabilization spot check (independent)
from buildintegrand import build_integrand_pf
def reduce_k_gm(k, want_keys):
    integrand_pf = build_integrand_pf(k)
    work = {}
    for key, R in integrand_pf.items():
        kk, RR = fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    out = {}
    for deg in range(k, 0, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            simple, higher, const = pf_split(R)
            if not higher[1] and not higher[2]: continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            if key in want_keys: out[key] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return out
g20_pf = pf_from_sympy((1+S2)/4 - S2/(2*(v**2+1)))
g11_pf = pf_from_sympy(-S2*(v**2+2*v-1)/(4*(v**2+1)))
def pf_nz(P):
    return bool(P[0]) or any(not e.is_zero() for e in P[1]) or any(not e.is_zero() for e in P[2]) or not P[3].is_zero()
t0=time.time()
out = reduce_k_gm(13, {(2,0),(1,1)})
print(f"k=13: Gm_20==g20: {not pf_nz(pf_sub(out[(2,0)], g20_pf))}, Gm_11==g11: {not pf_nz(pf_sub(out[(1,1)], g11_pf))} [{time.time()-t0:.0f}s]")
