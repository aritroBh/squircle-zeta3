"""WT5 step A2: control k=5 vs omega5.json with merged same-pole terms; deltas; cubic law."""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fastred
from fastred import (QSI, ZERO, AQ, BQ, pf_zero, pf_add, pf_scale, pf_sub)
from fractions import Fraction as Fr
from math import comb
t0 = time.time()

def pf_from_omega(terms):
    NAME = {"0": ZERO, "a": AQ, "-a": -AQ, "b": BQ, "-b": -BQ}
    P = pf_zero()
    for nm, (p, q) in terms:
        P[0][NAME[nm]] = P[0].get(NAME[nm], ZERO) + QSI(Fr(p), Fr(q))
    for r in list(P[0]):
        if P[0][r].is_zero(): del P[0][r]
    return P

def pf_to_omega(P):
    POLE_NAME = {ZERO: "0", AQ: "a", -AQ: "-a", BQ: "b", -BQ: "-b"}
    real, ipos, ineg, const = P
    assert not ipos and not ineg and const.is_zero()
    out = []
    for r, c in sorted(real.items(), key=lambda rc: ["0","a","-a","b","-b"].index(POLE_NAME[rc[0]])):
        assert c.is_real()
        out.append([POLE_NAME[r], [str(c.p0), str(c.p1)]])
    return out

def is_zero(diff):
    return (not diff[0]) and (not diff[1]) and (not diff[2]) and diff[3].is_zero()

kernels = json.load(open(os.path.join(HERE, "kernels_deg4.json")))
om5 = json.load(open(os.path.join(HERE, "omega5.json")))
def K(k, key): return pf_from_omega(kernels[str(k)][key])

# CONTROL merged
for key in ["4,0", "3,1", "2,2"]:
    diff = pf_sub(K(5, key), pf_from_omega(om5[key]))
    assert is_zero(diff), (key, pf_to_omega(diff))
print("CONTROL k=5 kernels == M5b omega5.json (4,0),(3,1),(2,2) after merging: EXACT MATCH", flush=True)

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

ok = True
for k in range(9, 15):
    for key in ["4,0", "3,1", "2,2"]:
        D1 = pf_from_omega(deltas[key]["D1"]); D2 = pf_from_omega(deltas[key]["D2"]); D3 = pf_from_omega(deltas[key]["D3"])
        pred = pf_add(pf_add(pf_add(K(5, key), pf_scale(D1, Fr(k-5))), pf_scale(D2, Fr(comb(k-5, 2)))), pf_scale(D3, Fr(comb(k-5, 3))))
        z = is_zero(pf_sub(pred, K(k, key)))
        ok = ok and z
        if not z: print(f"  LAW FAILS k={k} key={key}", flush=True)
print(f"CUBIC LAW k=9..14 all three keys exact: {ok}", flush=True)
assert ok
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_A2: k=5==omega5 (merged) EXACT; deltas_deg4.json; cubic law k=9..14 EXACT; t={time.time()-t0:.1f}s\n")
print("DONE", time.time()-t0, "s")
