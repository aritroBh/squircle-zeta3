import sys, importlib, pickle, os, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

# Extend record to k=14 with Gm for stabilization check of Gm20, Gm11
from fastred import pf_integrate_higher, pf_eval_at, pf_mul_DL, pf_split, pf_add, pf_scale
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

g20_pf   = pf_from_sympy((1+sp.sqrt(2))/4 - sp.sqrt(2)/(2*(v**2+1)))
# stabilization of Gm20, Gm11 to k=14
t0=time.time()
stab = True
gms = {}
for k in range(2, 15):
    out = reduce_k_gm(k, {(2,0),(1,1)})
    gms[k] = out
    if k >= 3:
        d20 = pf_sub(out[(2,0)], gms[k-1][(2,0)])
        z20 = (not d20[1] or all(e.is_zero() for e in d20[1])) and (not d20[2] or all(e.is_zero() for e in d20[2])) and d20[3].is_zero() and not d20[0]
        if not z20: stab=False; print(f"Gm20 NOT stable at k={k}")
    if k >= 4:
        d11 = pf_sub(out[(1,1)], gms[k-1][(1,1)])
        z11 = (not d11[1] or all(e.is_zero() for e in d11[1])) and (not d11[2] or all(e.is_zero() for e in d11[2])) and d11[3].is_zero() and not d11[0]
        if not z11: stab=False; print(f"Gm11 NOT stable at k={k}")
print(f"Gm20 stable for k=2..14 and Gm11 stable for k=3..14: {stab}  [{time.time()-t0:.0f}s]")
with open('progress.log','a') as f:
    f.write(f"\nStabilization Gm20 (k=2..14), Gm11 (k=3..14): {stab}\n")
    f.write("Conditional proof: (1-eps)Wcas_10 = eps^2 Rw2 + 0 eps^3  =>  S_10 = eps R10 + eps^2/(1-eps) Rw2\n")
