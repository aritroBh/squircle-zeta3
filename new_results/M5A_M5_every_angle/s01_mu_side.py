"""s01: formal coefficient P_{5,n} for symbolic n as Q(n)-combinations of Z-words Z_w(2n), per zeta monomial.
Uses mu_side() of g8_allN.py (copied unchanged from paper/ieee/proofs/scripts). Saves mu5_zsum.pkl.
Also tests against the independent numeric routine g6_common.mu_formal for n=1..8."""
import sys, time, pickle
sys.setrecursionlimit(20000)
from fractions import Fraction as Fr
import g8_allN as G
import g6_common as gc
t0 = time.time()
P = G.mu_side(5)
print(f"mu side k=5 built in {time.time()-t0:.1f}s")
for zm, E in sorted(P.items()):
    idx = sorted({a for w in E for a in w})
    print(f"  {G.zmon_name(zm):8s}: {len(E):4d} words, max depth {max((len(w) for w in E), default=0)}, indices {idx}")
ok = True
for nv in range(1, 9):
    tab = gc.mu_formal(nv, 5)[5]
    keys = set(tab) | set(P)
    for zm in keys:
        v = G.ev(P.get(zm, {}), nv)
        if v != tab.get(zm, Fr(0)):
            ok = False; print("MISMATCH", nv, zm, v, tab.get(zm))
print("mu side == g6_common.mu_formal (independent, numeric n), n=1..8:", ok)
# denominators
import sympy as sp
dens = set()
for zm, E in P.items():
    for w, c in E.items():
        d = sp.factor_list(c.denom.as_expr())[1]
        for f, m in d: dens.add((str(f), m))
print("denominator factors (factor, max mult):")
mx = {}
for f, m in dens: mx[f] = max(mx.get(f, 0), m)
print("  ", mx)
pickle.dump({zm: {w: (str(c.numer.as_expr()), str(c.denom.as_expr())) for w, c in E.items()} for zm, E in P.items()},
            open("mu5_zsum.pkl", "wb"))
