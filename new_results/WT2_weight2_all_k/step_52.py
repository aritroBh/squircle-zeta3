import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

# Recompute record with H for k up to 12 (fresh, to include higher degrees)
from fastred import pf_integrate_higher, pf_eval_at, pf_mul_DL, pf_split, pf_add, pf_scale, pf_from_poly
def reduce_k_instr(k, record=None):
    integrand_pf = build_integrand_pf(k)
    work = {}
    for key, R in integrand_pf.items():
        kk, RR = fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    final = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if record is not None:
                s0, h0, c0 = pf_split(R)
                record[(k, key)] = {'W': R, 'H': h0}
            if deg == 0:
                final[key] = R; continue
            simple, higher, const = pf_split(R)
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]: continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            if record is not None: record[(k,key)]['Gm'] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return final

record = {}
for k in range(1, 13):
    reduce_k_instr(k, record=record)
pickle.dump(record, open('record2.pkl','wb'))

# Express a pure-ipart PF as a real rational function in v (combine conjugates) via sympy
v = sp.Symbol('v', positive=True)
def ipart_to_real(P):
    return sp.together(pf_to_sympy([{}, P[1], P[2], ZERO]))

# stable H: print H^{(k)} for largest k and check stabilization window
print("Stable higher parts h_{ij}(v) (real form):")
for key in [(2,0),(1,1),(3,0),(2,1),(4,0),(3,1),(2,2),(5,0),(4,1),(3,2)]:
    d = sum(key)
    Hk = record[(12 if 12>=d else d, key)]['H'] if (12,key) in record else None
    # use k=12
    H = record[(12, key)]['H']
    # check stable: H same for k=11,12?
    H11 = record[(11, key)]['H']
    same = pf_sub(H, H11)
    stable = (not same[1] or all(e.is_zero() for e in same[1])) and (not same[2] or all(e.is_zero() for e in same[2]))
    print(f"{key}: stable(11=12):{stable}   h = {sp.simplify(ipart_to_real(H))}")
