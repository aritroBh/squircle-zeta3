"""WT3 step 1: degree-2 kernels R_20, R_11 for k=3..15 (fast reducer), control vs ibp_reduce at k=3,4."""
import sys, os, json, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
sys.path.insert(0, WORK)
sys.path.insert(0, os.path.join(BASE, 'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
from fastred import *
from buildintegrand import build_integrand_pf
import fastred

def pf_is_zero(P):
    return (not P[0]) and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

def qsi_str(q):
    assert q.is_real()
    return [str(q.p0), str(q.p1)]

def kernel_to_json(P):
    # only real simple poles expected; keep +-i simple-pole coeffs too if present
    assert P[3].is_zero(), "nonzero constant part"
    real = sorted(((r, c) for r, c in P[0].items()), key=lambda rc: (rc[0].p0, rc[0].p1))
    out = {"real": [[qsi_str(r), qsi_str(c)] for r, c in real]}
    if P[1] or P[2]:
        out["ipos"] = [[str(e.p0), str(e.p1), str(e.p2), str(e.p3)] for e in P[1]]
        out["ineg"] = [[str(e.p0), str(e.p1), str(e.p2), str(e.p3)] for e in P[2]]
    return out

results = {}
timings = {}
for k in range(3, 16):
    t0 = time.time()
    pf_int = build_integrand_pf(k)
    t1 = time.time()
    red = fastred.reduce_k_fast(k, integrand_pf=pf_int)
    t2 = time.time()
    deg2 = {key: red.get(key, pf_zero()) for key in [(2,0),(1,1)]}
    # sanity: no +-i poles, no const, real residues
    for key, P in deg2.items():
        pf_check_reality(P, ctx=f"k={k} {key}")
        assert not P[1] and not P[2], f"k={k} {key}: +-i poles survive: {P[1]} {P[2]}"
        assert P[3].is_zero()
    results[k] = {f"{i}{j}": kernel_to_json(P) for (i,j), P in deg2.items()}
    timings[k] = (round(t1-t0,2), round(t2-t1,2))
    print(f"k={k}: build {t1-t0:.1f}s reduce {t2-t1:.1f}s  poles20={sorted((str(r.p0),str(r.p1)) for r in deg2[(2,0)][0])}", flush=True)

with open(os.path.join(WORK, 'wt3_kernels.json'), 'w') as fh:
    json.dump({"kernels": {str(k): results[k] for k in results}, "timings": {str(k): timings[k] for k in timings}}, fh, indent=1)
print("saved wt3_kernels.json")

# ---------------- control: slow reference at k=3,4, ALL keys ----------------
import sympy as sp
from ibp_reduce import reduce_k
def pf_from_slow(R):
    return pf_from_sympy(R)
ok = True
for k in [3, 4]:
    t0 = time.time()
    slow = reduce_k(k, verbose=False)
    pf_int = build_integrand_pf(k)
    fast = fastred.reduce_k_fast(k, integrand_pf=pf_int)
    keys = set(slow) | set(fast)
    for key in sorted(keys, key=lambda kk: (-sum(kk), kk)):
        Ps = pf_from_slow(slow.get(key, 0))
        Pf = fast.get(key, pf_zero())
        same = pf_is_zero(pf_sub(Ps, Pf))
        if not same:
            ok = False
        print(f"  k={k} key={key}: {'MATCH' if same else 'MISMATCH'}", flush=True)
    print(f"k={k} control done in {time.time()-t0:.1f}s", flush=True)
print("CONTROL:", "PASS" if ok else "FAIL")
