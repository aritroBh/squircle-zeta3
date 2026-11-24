"""verify step: PREDICTIVE exact test of the polynomial-in-k law beyond the range used to find it.
Runs the fast reducer (read-only import of this folder's fastred/buildintegrand/genintegrand) for NEW k = 16, 17, ...,
extracts the weight-6 kernels (keys i+j = 5) and compares every residue (exact Q(sqrt2) rationals) with the
polynomial of degree <= 4 in k interpolated from k = 6..10 (exact Lagrange).  Also the same for weights 4,5,7,8
(degree d-1 polynomials in k, interpolated from the first d points k = d+1..2d).
usage: verify_extend_exact.py OUTDIR K [K ...]      (writes OUTDIR/ext_finals_K.json, resumable)
"""
import sys, os, json, time, pickle
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
from fractions import Fraction as Fr
import verify_lib as VL
from fastred import reduce_k_fast, pf_sub
from buildintegrand import build_integrand_pf

OUT = sys.argv[1]
KS = [int(x) for x in sys.argv[2:]]
os.makedirs(OUT, exist_ok=True)

def q2s(q):
    return [str(q.p0), str(q.p1), str(q.p2), str(q.p3)]

def pf2j(P):
    return {"real": [[q2s(r), q2s(c)] for r, c in sorted(P[0].items(), key=lambda rc: (rc[0].p0, rc[0].p1))],
            "ipos": [q2s(e) for e in P[1]], "ineg": [q2s(e) for e in P[2]], "const": q2s(P[3])}

def poly_through(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co

for k in KS:
    fn = os.path.join(OUT, f"ext_finals_{k}.json")
    if not os.path.exists(fn):
        t0 = time.time()
        ipf = build_integrand_pf(k)
        t1 = time.time()
        res = reduce_k_fast(k, integrand_pf=ipf)
        t2 = time.time()
        print(f"k={k}: build {t1-t0:.0f}s reduce {t2-t1:.0f}s", flush=True)
        data = {"k": k, "finals": {f"{i},{j}": pf2j(P) for (i, j), P in res.items()}}
        json.dump(data, open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
    new = VL.load_finals(fn)
    for d in (3, 4, 5, 6, 7):
        deg = d - 1
        pts_k = list(range(d + 1, d + 1 + deg + 1))          # first deg+1 points of the regime k >= d+1
        fins = {kk: VL.load_finals(os.path.join(W, f"finals_{kk}.json")) for kk in pts_k}
        keys = sorted({key for kk in pts_k for key in fins[kk] if sum(key) == d} | {key for key in new if sum(key) == d})
        nbad = 0; ncomp = 0
        for key in keys:
            labs = set(new.get(key, {}))
            for kk in pts_k:
                labs |= set(fins[kk].get(key, {}))
            for lab in sorted(labs):
                for comp in (0, 1):
                    ys = [(Fr(kk), fins[kk].get(key, {}).get(lab, (Fr(0), Fr(0)))[comp]) for kk in pts_k]
                    co = poly_through(ys)
                    pred = sum(co[j] * Fr(k) ** j for j in range(len(co)))
                    act = new.get(key, {}).get(lab, (Fr(0), Fr(0)))[comp]
                    ncomp += 1
                    if pred != act:
                        nbad += 1
                        print(f"   k={k} d={d} key {key} pole {lab} comp {comp}: predicted {pred} actual {act}")
        print(f"k={k}: weight {d+1}: {ncomp} residue components compared with the degree-{deg} polynomial from k={pts_k[0]}..{pts_k[-1]}: "
              f"mismatches = {nbad}", flush=True)
