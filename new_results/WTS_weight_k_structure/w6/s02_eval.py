"""w6 step 2: numerical value of the weight-(d+1) part W_{d+1}(k) from the exact final pieces.

usage: s02_eval.py K D [D ...] [--dps N]
   W_{d+1}(k) = sum_{i+j=d, i>=j} int_b^a K_ij(v) log^i t log^j(1-t) dv
using finals_K.json (exact, from s01_reduce.py) and tanh-sinh (mpmath) split at v = 1.
To keep the endpoint logarithms exact we integrate in the distance to the endpoint:
  left half  v = b + s, s in (0, 1-b):  t = s (v+a)/(4v),  log t, log(1-t) = log1p(-t)
  right half v = a - u, u in (0, a-1):  1-t = u (v+b)/(4v), log(1-t), log t = log1p(-(1-t))
(no cancellation near v=b or v=a; a simple pole at b (resp. a) is 1/s (resp. -1/u)).
Writes vals/k{K}_d{D}.json (dps 70) or vals_<dps>/... (resumable: skips existing files).
"""
import sys, os, json, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
from fractions import Fraction as Fr
import mpmath as mp

args = [a for a in sys.argv[1:]]
DPS = 70
if "--dps" in args:
    n = args.index("--dps"); DPS = int(args[n + 1]); del args[n:n + 2]
K = int(args[0]); DS = [int(x) for x in args[1:]]
VALS = "vals" if DPS == 70 else f"vals_{DPS}"
os.makedirs(os.path.join(W, VALS), exist_ok=True)

# precision is set AFTER all imports
mp.mp.dps = DPS
print(f"mp.mp.dps = {mp.mp.dps}", flush=True)
r2 = mp.sqrt(2)
aN = 1 + r2
bN = r2 - 1

def frac(s):
    f = Fr(s)
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)

def real_of(q):
    """QSI strings [p0,p1,p2,p3] -> mpf; asserts real"""
    assert Fr(q[2]) == 0 and Fr(q[3]) == 0, q
    return frac(q[0]) + frac(q[1]) * r2

data = json.load(open(os.path.join(W, f"finals_{K}.json")))
finals = data["finals"]

def terms_of(key):
    P = finals.get(key)
    if P is None:
        return None
    assert not P["ipos"] and not P["ineg"], f"{key}: +-i poles survive"
    assert all(Fr(x) == 0 for x in P["const"]), f"{key}: nonzero constant"
    out = []
    for r, c in P["real"]:
        rr = (Fr(r[0]), Fr(r[1]))
        tag = "0" if rr == (0, 0) else "a" if rr == (1, 1) else "b" if rr == (-1, 1) else "other"
        out.append((tag, real_of(r), real_of(c)))
    return out

def piece(i, j, ts):
    # convergence check on the endpoint poles
    for tag, r, c in ts:
        if tag == "b":
            assert j >= 1, f"pole at b against j=0 in key {(i,j)}: divergent"
        if tag == "a":
            assert i >= 1, f"pole at a against i=0 in key {(i,j)}: divergent"
    def K_left(s):
        v = bN + s
        tot = mp.mpf(0)
        for tag, r, c in ts:
            tot += c / (s if tag == "b" else (v - r))
        return v, tot
    def K_right(u):
        v = aN - u
        tot = mp.mpf(0)
        for tag, r, c in ts:
            tot += c / (-u if tag == "a" else (v - r))
        return v, tot
    def fl(s):
        v, Kv = K_left(s)
        t = s * (v + aN) / (4 * v)
        return Kv * mp.log(t) ** i * mp.log1p(-t) ** j
    def fr(u):
        v, Kv = K_right(u)
        w = u * (v + bN) / (4 * v)
        return Kv * mp.log1p(-w) ** i * mp.log(w) ** j
    v1, e1 = mp.quad(fl, [0, 1 - bN], error=True, maxdegree=12)
    v2, e2 = mp.quad(fr, [0, aN - 1], error=True, maxdegree=12)
    return v1 + v2, max(e1, e2)

for d in DS:
    fn = os.path.join(W, VALS, f"k{K}_d{d}.json")
    if os.path.exists(fn):
        print(f"k={K} d={d}: exists, skipping", flush=True)
        continue
    t0 = time.time()
    total = mp.mpf(0)
    emax = mp.mpf(0)
    per = {}
    for i in range(d, (d - 1) // 2, -1):
        j = d - i
        ts = terms_of(f"{i},{j}")
        if ts is None:
            continue
        val, err = piece(i, j, ts)
        per[f"{i},{j}"] = {"value": mp.nstr(val, DPS), "err": mp.nstr(err, 5),
                           "poles": [t[0] for t in ts]}
        total += val
        emax = max(emax, err)
        print(f"  k={K} d={d} key=({i},{j}) val={mp.nstr(val, 30)} err={mp.nstr(err, 3)} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    rec = {"k": K, "d": d, "weight": d + 1, "dps": mp.mp.dps, "value": mp.nstr(total, DPS),
           "err_est": mp.nstr(emax, 5), "pieces": per, "seconds": time.time() - t0}
    tmp = fn + ".tmp"
    json.dump(rec, open(tmp, "w"), indent=1)
    os.replace(tmp, fn)
    print(f"k={K} W_{d+1} = {mp.nstr(total, 50)}  err_est {mp.nstr(emax, 3)}  dps {mp.mp.dps}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
