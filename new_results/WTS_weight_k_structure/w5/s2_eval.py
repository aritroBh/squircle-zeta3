# Step 2: numerical evaluation of the weight-(d+1) final pieces,
#   W_{d+1}(k) = sum_{i+j=d} int_b^a K_ij^{(k)}(v) log^i t log^j(1-t) dv,
# from the exact final pieces final_<k>.pkl (written by s1_reduce.py).
# usage: s2_eval.py <deg d> <dps> <k> [<k> ...]      (tanh-sinh, split at v=1)
# Each result is saved to results/w<d+1>_k<k>_dps<dps>.json (resumable: existing files are skipped).
#
# Endpoint treatment: on [b,1] use v = b+s (s = v-b exact), on [1,a] use v = a-s (s = a-v exact),
# so the factorised forms  t = (v-b)(v+a)/(4v),  1-t = (a-v)(v+b)/(4v)  and the pole factors
# v-b, v-a are evaluated without cancellation near the singular endpoints.
import sys, os, time, json, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import fastred
from fastred import QSI, ZERO, AQ, BQ
from fractions import Fraction as Fr

d = int(sys.argv[1]); dps = int(sys.argv[2]); ks = [int(x) for x in sys.argv[3:]]
mp.mp.dps = dps            # AFTER all imports (mp.mp.dps, not mp.dps)
print("mp.mp.dps =", mp.mp.dps, flush=True)

def fr(x): return mp.mpf(x.numerator) / x.denominator
def q2mpf(q):
    assert q.p2 == 0 and q.p3 == 0, q
    return fr(q.p0) + fr(q.p1) * mp.sqrt(2)

def evaluate(k):
    final = pickle.load(open(f"final_{k}.pkl", "rb"))
    keys = sorted(kk for kk in final if sum(kk) == d)
    sq2 = mp.sqrt(2); a = 1 + sq2; b = sq2 - 1
    # pole labels
    labels = {}
    for name, q in (("0", ZERO), ("a", AQ), ("-a", -AQ), ("b", BQ), ("-b", -BQ)):
        labels[(q.p0, q.p1)] = name
    terms = []   # (i, j, pole_label, coeff)
    for key in keys:
        P = final[key]
        assert not P[1] and not P[2] and P[3].is_zero(), ("non-real-letter content", key)
        for r, c in P[0].items():
            lab = labels[(r.p0, r.p1)]
            if lab == "b":
                assert key[1] >= 1, ("pole at b with j=0 (divergent)", key)
            terms.append((key[0], key[1], lab, q2mpf(c)))
    # collect by (i,j)
    ij = sorted(set((i, j) for i, j, _, _ in terms))
    def make(side):
        def f(s):
            if side == 1:           # v = b + s
                v = b + s
                t = s * (v + a) / (4 * v)
                omt = (2 - s) * (v + b) / (4 * v)
                l1 = mp.log(t); l2 = mp.log1p(-t)
                dv = {"0": v, "a": s - 2, "-a": v + a, "b": s, "-b": v + b}
            else:                   # v = a - s
                v = a - s
                omt = s * (v + b) / (4 * v)
                t = (2 - s) * (v + a) / (4 * v)
                l2 = mp.log(omt); l1 = mp.log1p(-omt)
                dv = {"0": v, "a": -s, "-a": v + a, "b": 2 - s, "-b": v + b}
            tot = mp.mpf(0)
            for i, j, lab, c in terms:
                tot += c / dv[lab] * l1**i * l2**j
            return tot
        return f
    t0 = time.time()
    I1, e1 = mp.quad(make(1), [0, 1 - b], error=True, maxdegree=12)
    I2, e2 = mp.quad(make(2), [0, a - 1], error=True, maxdegree=12)
    return I1 + I2, e1 + e2, time.time() - t0, keys

for k in ks:
    fn = f"results/w{d+1}_k{k}_dps{dps}.json"
    if os.path.exists(fn):
        print(f"k={k}: {fn} exists, skipping", flush=True); continue
    val, err, dt, keys = evaluate(k)
    rec = dict(weight=d + 1, k=k, dps=dps, value=mp.nstr(val, dps), quad_err_estimate=mp.nstr(err, 5),
               seconds=dt, keys=[list(x) for x in keys])
    json.dump(rec, open(fn, "w"), indent=1)
    print(f"k={k} W_{d+1}={mp.nstr(val, 40)} err~{mp.nstr(err, 3)} ({dt:.1f}s) dps={mp.mp.dps}", flush=True)
