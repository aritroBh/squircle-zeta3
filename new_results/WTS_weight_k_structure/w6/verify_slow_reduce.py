"""verify step: independent exact reduction of the M_k integrand by the SLOW sympy reducer
(paper/squircle_note/scripts/ibp_reduce.reduce_k, whose integrand comes from mk_integrands.v_integrand:
symbolic d^k/de^k of  A e^{e l1} + B e^{e l2} - sqrt(A^2 e^{2 e l1} + B^2 e^{2 e l2}), i.e. NOT from the
closed form genintegrand/buildintegrand used by the fast reducer).

usage: verify_slow_reduce.py K OUTDIR
Writes OUTDIR/slow_finals_K.json : {key: {pole_label: [p0,p1]}} with exact Fractions as strings,
pole labels in {0,a,-a,b,-b,i,-i}; residue extraction is cross-checked by exact re-summation
(sum of residue/(v-r) minus the sympy expression must simplify to 0).
"""
import sys, os, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v, S2

K = int(sys.argv[1]); OUT = sys.argv[2]
os.makedirs(OUT, exist_ok=True)
a_, b_ = 1 + S2, S2 - 1
POLES = {"0": sp.Integer(0), "a": a_, "-a": -a_, "b": b_, "-b": -b_}

def q2(x):
    """x in Q(sqrt2) -> [p0, p1] strings"""
    x = sp.expand(sp.radsimp(x))
    p0 = x.subs(S2, 0)
    p1 = sp.expand((x - p0) / S2)
    assert not p0.free_symbols and not p1.free_symbols, x
    return [str(sp.Rational(p0)), str(sp.Rational(p1))]

t0 = time.time()
fin = reduce_k(K, verbose=False)
print(f"k={K}: slow reduce done in {time.time()-t0:.1f}s; keys {sorted(fin, key=lambda kk: (-sum(kk), kk))}", flush=True)
res = {}
for key, expr in fin.items():
    R = sp.together(expr)
    rec = {}
    tot = 0
    for lab, r in POLES.items():
        c = sp.cancel(sp.together(R * (v - r)))
        c = sp.radsimp(sp.simplify(c.subs(v, r)))
        if c != 0:
            rec[lab] = q2(c)
            tot += c / (v - r)
    # check: nothing left (no other poles, no polynomial part, no +-i poles)
    left = sp.simplify(sp.together(R - tot))
    rec["_leftover_is_zero"] = (left == 0)
    if left != 0:
        rec["_leftover"] = str(left)
    res[f"{key[0]},{key[1]}"] = rec
    print(f"  key {key}: poles {[l for l in rec if not l.startswith('_')]} leftover0={left == 0}", flush=True)
json.dump({"k": K, "slow": res, "seconds": time.time() - t0}, open(os.path.join(OUT, f"slow_finals_{K}.json"), "w"), indent=1)
print("saved", flush=True)
