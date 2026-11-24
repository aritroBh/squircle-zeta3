import sys, importlib, pickle, os, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

src='<project>/research/08_theory/WT2_weight2_all_k'
# reuse cached binteg pkls from WT2 to save time (they are the exact integrands, already verified)
def get_integrand(k):
    fn = f"{src}/binteg_{k}.pkl"
    if os.path.exists(fn):
        return pickle.load(open(fn,'rb'))
    return build_integrand_pf(k)

def reduce_k_instr(k, record=None):
    integrand_pf = get_integrand(k)
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
            assert const.is_zero()
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]:
                if record is not None: record[(k,key)]['Gm'] = pf_zero()
                continue
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
t0=time.time()
for k in range(1, 13):
    reduce_k_instr(k, record=record)
print(f"computed records k<=12 in {time.time()-t0:.1f}s")
pickle.dump(record, open('record_W2b.pkl','wb'))

# verify stabilization of Gm_20 (k=2..12) and Gm_11 (k=3..12)
def pf_eq(P, Q):
    d = pf_sub(P, Q)
    return (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
g20_ref = record[(2,(2,0))]['Gm']
stab20 = all(pf_eq(record[(k,(2,0))]['Gm'], g20_ref) for k in range(2,13))
g11_ref = record[(3,(1,1))]['Gm']
stab11 = all(pf_eq(record[(k,(1,1))]['Gm'], g11_ref) for k in range(3,13))
print("Gm20 stable k=2..12:", stab20, " Gm11 stable k=3..12:", stab11)
