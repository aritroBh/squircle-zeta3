base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step A: exact degree-4 kernels K^(k) for k=5..14 via fastred; omega-style JSON;
D1 = K6-K5, D2 = K7-2K6+K5, D3 = K8-3K7+3K6-K5; cubic law check k=9..14; control vs M5b omega5.json."""
import sys, os, json, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fastred
from fastred import (QSI, ZERO, AQ, BQ, reduce_k_fast, pf_zero, pf_add, pf_scale, pf_sub)
from fractions import Fraction as Fr

t0 = time.time()
POLE_NAME = {ZERO: "0", AQ: "a", -AQ: "-a", BQ: "b", -BQ: "-b"}
KEYS = [(4, 0), (3, 1), (2, 2)]

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

path = os.path.join(HERE, "kernels_deg4.json")
kernels = json.load(open(path)) if os.path.exists(path) else {}
for k in range(5, 15):
    if str(k) in kernels:
        print(f"k={k}: loaded from checkpoint", flush=True)
        continue
    ipf = pickle.load(open(os.path.join(HERE, f"binteg_{k}.pkl"), "rb"))
    res = reduce_k_fast(k, integrand_pf=ipf)
    rec = {}
    for key in KEYS:
        rec[",".join(map(str, key))] = pf_to_omega(res[key])
    kernels[str(k)] = rec
    json.dump(kernels, open(path, "w"), indent=1)
    print(f"k={k}: (4,0) terms={len(rec['4,0'])} (3,1) terms={len(rec['3,1'])} (2,2) terms={len(rec['2,2'])}  t={time.time()-t0:.1f}s", flush=True)

# CONTROL: k=5 vs M5b omega5.json keys (4,0),(3,1),(2,2)
om5 = json.load(open(os.path.join(HERE, "omega5.json")))
for key in ["4,0", "3,1", "2,2"]:
    mine = sorted((t[0], t[1][0], t[1][1]) for t in kernels["5"][key])
    theirs = sorted((t[0], t[1][0], t[1][1]) for t in om5[key])
    assert mine == theirs, (key, mine, theirs)
print("CONTROL k=5 kernels == M5b omega5.json (4,0),(3,1),(2,2): EXACT MATCH", flush=True)

def K(k, key): return pf_from_omega(kernels[str(k)][key])

# D1,D2,D3 as PF and omega
deltas = {}
for key in ["4,0", "3,1", "2,2"]:
    D1 = pf_sub(K(6, key), K(5, key))
    D2 = pf_add(pf_sub(K(7, key), pf_scale(K(6, key), Fr(2))), K(5, key))
    D3 = pf_add(pf_add(pf_sub(K(8, key), pf_scale(K(7, key), Fr(3))), pf_scale(K(6, key), Fr(3))), pf_scale(K(5, key), Fr(-1)))
    deltas[key] = {"D1": pf_to_omega(D1), "D2": pf_to_omega(D2), "D3": pf_to_omega(D3)}
    print(f"{key}: D1={deltas[key]['D1']}", flush=True)
    print(f"{key}: D2={deltas[key]['D2']}", flush=True)
    print(f"{key}: D3={deltas[key]['D3']}", flush=True)
json.dump(deltas, open(os.path.join(HERE, "deltas_deg4.json"), "w"), indent=1)

# CUBIC LAW: K^(k) == K5 + (k-5) D1 + binom(k-5,2) D2 + binom(k-5,3) D3 for k=9..14 (exact)
from math import comb
ok = True
for k in range(9, 15):
    for key in ["4,0", "3,1", "2,2"]:
        D1 = pf_from_omega(deltas[key]["D1"]); D2 = pf_from_omega(deltas[key]["D2"]); D3 = pf_from_omega(deltas[key]["D3"])
        pred = pf_add(pf_add(pf_add(K(5, key), pf_scale(D1, Fr(k-5))), pf_scale(D2, Fr(comb(k-5, 2)))), pf_scale(D3, Fr(comb(k-5, 3))))
        diff = pf_sub(pred, K(k, key))
        z = (not diff[0]) and (not diff[1]) and (not diff[2]) and diff[3].is_zero()
        ok = ok and z
        if not z: print(f"  LAW FAILS k={k} key={key}", flush=True)
print(f"CUBIC LAW k=9..14 all three keys exact: {ok}", flush=True)
assert ok
# k=4 OFF the cubic? compute k=4 kernels too for the offset record (binteg_4 from WT4 not copied; skip)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_A kernels: k=5..14 reduced; k=5==omega5 EXACT; cubic law k=9..14 EXACT; t={time.time()-t0:.1f}s\\n")
print("DONE", time.time()-t0, "s")
'''
open(W+'wt5_A_kernels.py','w').write(script)
print('written', len(script))
