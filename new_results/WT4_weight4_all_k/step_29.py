import os, textwrap
work='<project>/research/08_theory/WT4_weight4_all_k'
script = r'''"""WT4 step A: exact degree-3 kernels K^(k) for k=3..14 via fastred; omega-style JSON;
Delta1 = K5-K4, Delta2 = K6-2K5+K4; quadratic law check for k=7..14; control vs M4a omega4.json."""
import sys, os, json, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fastred
from fastred import (QSI, ZERO, AQ, BQ, reduce_k_fast, pf_zero, pf_add, pf_scale, pf_sub)
from fractions import Fraction as Fr

t0 = time.time()
POLE_NAME = {ZERO: "0", AQ: "a", -AQ: "-a", BQ: "b", -BQ: "-b"}

def pf_to_omega(P):
    real, ipos, ineg, const = P
    assert not ipos and not ineg, ("+-i poles survive", ipos, ineg)
    assert const.is_zero(), ("nonzero const", const)
    out = []
    for r, c in sorted(real.items(), key=lambda rc: ["0","a","-a","b","-b"].index(POLE_NAME[rc[0]])):
        assert c.is_real(), ("nonreal coefficient", r, c)
        out.append([POLE_NAME[r], [str(c.p0), str(c.p1)]])
    return out

def pf_from_omega(terms):
    NAME = {"0": ZERO, "a": AQ, "-a": -AQ, "b": BQ, "-b": -BQ}
    P = pf_zero()
    for nm, (p, q) in terms:
        P[0][NAME[nm]] = P[0].get(NAME[nm], ZERO) + QSI(Fr(p), Fr(q))
    return P

kernels = {}
for k in range(3, 15):
    ipf = pickle.load(open(os.path.join(HERE, f"binteg_{k}.pkl"), "rb"))
    res = reduce_k_fast(k, integrand_pf=ipf)
    rec = {}
    for key in [(2, 1), (3, 0)]:
        rec[",".join(map(str, key))] = pf_to_omega(res[key])
    kernels[k] = rec
    print(f"k={k}: (2,1) terms={len(rec['2,1'])} (3,0) terms={len(rec['3,0'])}  t={time.time()-t0:.1f}s", flush=True)

json.dump(kernels, open(os.path.join(HERE, "kernels_deg3.json"), "w"), indent=1)

# CONTROL: k=4 vs M4a omega4.json
M4A = os.path.join(HERE, "..", "M4a_decomposition")
om4 = json.load(open(os.path.join(M4A, "omega4.json")))
for key in ["2,1", "3,0"]:
    mine = sorted((t[0], t[1][0], t[1][1]) for t in kernels[4][key])
    theirs = sorted((t[0], t[1][0], t[1][1]) for t in om4[key])
    assert mine == theirs, (key, mine, theirs)
print("CONTROL k=4 kernels == M4a omega4.json (2,1),(3,0): EXACT MATCH", flush=True)

def K(k, key): return pf_from_omega(kernels[k][key])

# Delta1, Delta2 as PF and omega
deltas = {}
for key in ["2,1", "3,0"]:
    D1 = pf_sub(K(5, key), K(4, key))
    D2 = pf_add(pf_sub(K(6, key), pf_scale(K(5, key), Fr(2))), K(4, key))
    deltas[key] = {"D1": pf_to_omega(D1), "D2": pf_to_omega(D2)}
    print(f"{key}: D1={deltas[key]['D1']}  D2={deltas[key]['D2']}", flush=True)
json.dump(deltas, open(os.path.join(HERE, "deltas_deg3.json"), "w"), indent=1)

# QUADRATIC LAW: K^(k) == K4 + (k-4) D1 + binom(k-4,2) D2 for k=7..14 (exact)
from math import comb
ok = True
for k in range(7, 15):
    for key in ["2,1", "3,0"]:
        D1 = pf_from_omega(deltas[key]["D1"]); D2 = pf_from_omega(deltas[key]["D2"])
        pred = pf_add(pf_add(K(4, key), pf_scale(D1, Fr(k-4))), pf_scale(D2, Fr(comb(k-4, 2))))
        diff = pf_sub(pred, K(k, key))
        z = (not diff[0]) and (not diff[1]) and (not diff[2]) and diff[3].is_zero()
        ok = ok and z
        if not z: print(f"  LAW FAILS k={k} key={key}", flush=True)
print(f"QUADRATIC LAW k=7..14 both keys exact: {ok}", flush=True)
assert ok
# also k=5,6 by construction trivially; report k=3 offset later
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_A kernels: k=3..14 reduced; k=4==omega4 EXACT; quadratic law k=7..14 EXACT; t={time.time()-t0:.1f}s\n")
print("DONE", time.time()-t0, "s")
'''
p=os.path.join(work,'wt4_A_kernels.py')
open(p,'w').write(script)
print('written', p)
