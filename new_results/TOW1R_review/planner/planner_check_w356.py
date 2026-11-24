# Planner check at a NEW k (reducer data so far: k <= 21): one run of the paper's fast reducer, then
#  (a) weight 3: kernels R20, R11 vs TOW1/WT3 closed forms;
#  (b) weight 5: degree-4 kernels vs the M1 session's exact cubic residue polynomials (deg4_residue_polynomials.json);
#  (c) weight 6: degree-5 kernels vs its exact quartic residue polynomials (deg5_residue_polynomials.json);
#  (d) the weight-3 slope constant found by PSLQ, not assumed.
import sys, time, json, pickle
sys.dont_write_bytecode = True
R = "<project>/research/08_theory"
sys.path.insert(0, R + "/WTS_weight_k_structure/theory/")
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
from fastred import *
from buildintegrand import build_integrand_pf

k = int(sys.argv[1]) if len(sys.argv) > 1 else 23
t0 = time.time()
fin = reduce_k_fast(k, build_integrand_pf(k))
pickle.dump(fin, open(f"{R}/TOW1R_review/planner/fin_k{k}.pkl", "wb"))
print(f"reduced k={k} in {time.time()-t0:.0f}s; keys {sorted(fin)}", flush=True)

def got(P):
    d = {(r.p0, r.p1): (c.p0, c.p1) for r, c in P[0].items() if c.p0 or c.p1}
    imag = any(c.p2 or c.p3 for c in P[0].values())
    return d, imag

# (a) weight 3
r20 = {(0, 0): (Fr(1, 4), Fr(-1, 4)), (1, 1): (Fr(1, 4), 0), (-1, -1): (Fr(2*k-7, 4), 0), (1, -1): (Fr(-(k-2), 2), 0)}
r11 = {(0, 0): (Fr(2*k-5, 4), Fr(1, 4)), (-1, -1): (Fr(-(2*k-5), 2), 0)}
for key, exp in [((2, 0), r20), ((1, 1), r11)]:
    g, im = got(fin[key])
    print(f"(a) k={k} key {key}: match={g == exp} imag={im} higher={fin[key][1:]}", flush=True)

# (b), (c) weights 5 and 6
ROOT = {"0": (0, 0), "a": (1, 1), "-a": (-1, -1), "-b": (1, -1), "b": (-1, 1)}
K = sp.Symbol("k")
for tag, fn in [("(b) w5", "deg4_residue_polynomials.json"), ("(c) w6", "deg5_residue_polynomials.json")]:
    J = json.load(open(f"{R}/TOWN1_tower_large_n/claude_checks/{fn}"))
    for skey, roots in J.items():
        key = (int(skey[0]), int(skey[1]))
        exp = {}
        for rname, comp in roots.items():
            c0 = sp.Rational(sp.sympify(comp.get("rat", "0")).subs(K, k))
            c1 = sp.Rational(sp.sympify(comp.get("sqrt2", "0")).subs(K, k))
            if c0 or c1:
                exp[ROOT[rname]] = (Fr(int(c0.p), int(c0.q)), Fr(int(c1.p), int(c1.q)))
        g, im = got(fin[key])
        ok = g == exp
        print(f"{tag} k={k} key {key}: match={ok} imag={im}", flush=True)
        if not ok:
            print("   reducer:", g, "\n   json:   ", exp, flush=True)

# (d) slope constant
mp.mp.dps = 110
s2 = mp.sqrt(2); a = 1 + s2; b = s2 - 1
t = lambda v: (v*v + 2*v - 1)/(4*v)
om = lambda v: (1 + 2*v - v*v)/(4*v)
S20 = lambda v: 1/(2*(v + a)) - 1/(2*(v + b))
S11 = lambda v: 1/(2*v) - 1/(v + a)
slope = mp.quad(lambda v: S20(v)*mp.log(t(v))**2 + S11(v)*mp.log(t(v))*mp.log(om(v)), [b, 1, a])
rel = mp.pslq([slope, mp.zeta(3), mp.pi**2*mp.log(2), mp.log(2)**3], maxcoeff=10**6, maxsteps=10**6)
print("(d) PSLQ [slope, z3, pi^2 log2, log^3 2]:", rel, "(claim -> [24, 3, 1, 2])", flush=True)
