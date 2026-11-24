step=r'''import sys, importlib, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf
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
S2=sp.sqrt(2)
g20_pf = pf_from_sympy((1+S2)/4 - S2/(2*(v**2+1)))
g11_pf = pf_from_sympy(-S2*(v**2+2*v-1)/(4*(v**2+1)))
ng20 = pf_neg(g20_pf)
def pf_nz(P):
    return bool(P[0]) or any(not e.is_zero() for e in P[1]) or any(not e.is_zero() for e in P[2]) or not P[3].is_zero()
stab=True
for k in [13,14]:
    t0=time.time()
    out = reduce_k_gm(k, {(2,0),(1,1)})
    ok20 = not pf_nz(pf_sub(out[(2,0)], g20_pf))
    ok11 = not pf_nz(pf_sub(out[(1,1)], g11_pf))
    print(f"k={k}: Gm_20 == g20: {ok20},  Gm_11 == g11: {ok11}   [{time.time()-t0:.0f}s]")
    stab &= ok20 and ok11
print("Stabilization Lemma verified exactly for k=13,14 (extending k<=12):", stab)
with open('progress.log','a') as f:
    f.write(f"WT2b final: stabilization re-verified k=13,14 in this dir: {stab}\n")
'''
open('step_43.py','w').write(step)
import subprocess
r=subprocess.run(['<project>/.venv/bin/python','step_43.py'],capture_output=True,text=True,timeout=600)
print(r.stdout); print(r.stderr[-500:] if r.stderr else '')
