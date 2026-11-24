"""w6 step 1: exact fast reduction of the M_k integrand for one k.

usage: s01_reduce.py K
Reads binteg_K.pkl (copied from WT2_weight2_all_k), also rebuilds the integrand from the
closed form (build_integrand_pf) and checks the pickle equals it exactly (control), runs
reduce_k_fast, and writes finals_K.json: every reduced key (i,j) -> PF with exact Fraction
strings.  Resumable: does nothing if finals_K.json exists.
"""
import sys, os, json, pickle, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
from fastred import *            # QSI, pf_*, reduce_k_fast
from buildintegrand import build_integrand_pf

k = int(sys.argv[1])
out = os.path.join(W, f"finals_{k}.json")
if os.path.exists(out):
    print(f"k={k}: finals_{k}.json exists, skipping")
    sys.exit(0)

def q2s(q):
    return [str(q.p0), str(q.p1), str(q.p2), str(q.p3)]

def pf2j(P):
    return {"real": [[q2s(r), q2s(c)] for r, c in sorted(P[0].items(), key=lambda rc: (rc[0].p0, rc[0].p1))],
            "ipos": [q2s(e) for e in P[1]],
            "ineg": [q2s(e) for e in P[2]],
            "const": q2s(P[3])}

def pf_is_zero(P):
    return (not P[0]) and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

t0 = time.time()
ipf = pickle.load(open(os.path.join(W, f"binteg_{k}.pkl"), "rb"))
t1 = time.time()
ipf2 = build_integrand_pf(k)
t2 = time.time()
same = (set(ipf) == set(ipf2)) and all(pf_is_zero(pf_sub(ipf[key], ipf2[key])) for key in ipf)
print(f"k={k}: pickle == rebuilt closed-form integrand: {same}  (load {t1-t0:.1f}s, rebuild {t2-t1:.1f}s)", flush=True)
assert same
res = reduce_k_fast(k, integrand_pf=ipf)
t3 = time.time()
print(f"k={k}: reduce {t3-t2:.1f}s; keys: {sorted(res, key=lambda kk: (-sum(kk), kk))}", flush=True)
data = {"k": k, "same_as_rebuilt": same, "reduce_seconds": t3 - t2,
        "finals": {f"{i},{j}": pf2j(P) for (i, j), P in res.items()}}
tmp = out + ".tmp"
json.dump(data, open(tmp, "w"), indent=1)
os.replace(tmp, out)
print(f"k={k}: saved {out}", flush=True)
