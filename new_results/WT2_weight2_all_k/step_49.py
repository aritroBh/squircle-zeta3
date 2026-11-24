import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

def reduce_k_instr(k, integrand_pf=None, record=None):
    if integrand_pf is None:
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
                simple0, higher0, const0 = pf_split(R)
                rec = {'W': R, 'H': higher0}
                record[(k, key)] = rec
            if deg == 0:
                final[key] = R
                continue
            simple, higher, const = pf_split(R)
            assert const.is_zero()
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]:
                if record is not None: rec['Gm'] = pf_zero(); rec['C'] = ZERO
                continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            if record is not None:
                rec['G'] = G; rec['C'] = C; rec['Gm'] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return final

record = {}
N = 10
for k in range(1, N+1):
    reduce_k_instr(k, record=record)

# Series for Gm_{20} and Gm_{11}: look at pole structure across k
def show_series(keyname, key, upto=10):
    print(f"===== Gm series for key {keyname}")
    for k in range(2, upto+1):
        rec = record.get((k, key))
        if rec is None or 'Gm' not in rec: continue
        Gm = rec['Gm']
        parts = []
        for r, c in sorted(Gm[0].items(), key=lambda rc: float(rc[0].sympy())):
            parts.append(f"[{c.sympy()}]@({r.sympy()})")
        for jj, e in enumerate(Gm[1], start=1):
            if not e.is_zero(): parts.append(f"[{e.sympy()}]/(v-i)^{jj}")
        if not Gm[3].is_zero(): parts.append(f"const {Gm[3].sympy()}")
        print(f"  k={k}: " + " ".join(parts))

show_series("(2,0)", (2,0))
show_series("(1,1)", (1,1))
pickle.dump(record, open('record.pkl','wb'))
