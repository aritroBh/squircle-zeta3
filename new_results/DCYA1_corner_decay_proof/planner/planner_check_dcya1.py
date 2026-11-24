# Planner check of DCYA1's M_k(2pi/3) table: k <= 8 against the exact flat-angle series (ZLAD1's mu table, n <= 150, c^2 = 1/4),
# and the leading decay law at k = 40..80 (ratio to -(sqrt(2a sin a)/N) e^{-sqrt(2aN)} cos sqrt(2aN)).
import json, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
tab = json.load(open("<project>/research/08_theory/DCYA1_corner_decay_proof/mk_2pi3_k80.json"))
mu = json.load(open("<project>/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N150_K8.json"))
z = {f"z{w}": mp.zeta(w) for w in range(2, 9)}
def mono_val(m):
    if m == "1": return mp.mpf(1)
    v = mp.mpf(1)
    for part in m.replace("z", " z").split():
        if "^" in part: b, e = part.split("^"); v *= z[b] ** int(e)
        else: v *= z[part]
    return v
c2 = mp.mpf(1) / 4
val = lambda x: mp.mpf(x) if not isinstance(x, (list, dict)) else mp.mpf(x.get("value", x.get("M", "nan")) if isinstance(x, dict) else x[-1])
print("table entry type:", type(tab[0]).__name__, str(tab[0])[:80])
for k in range(0, 9):
    M = mp.fsum(mp.mpf(Fr(q).numerator) / Fr(q).denominator * mono_val(m) * c2 ** int(n) for n, terms in mu[str(k)].items() for m, q in terms.items())
    print(f"k={k}: flat series {mp.nstr(M, 25)} | DCYA1 table {mp.nstr(val(tab[k]), 25)} | diff {mp.nstr(abs(M - val(tab[k])), 3)}")
a = 2 * mp.pi / 3
for k in (40, 60, 80):
    N = k + 1; th = mp.sqrt(2 * a * N)
    lead = -(mp.sqrt(2 * a * mp.sin(a)) / N) * mp.exp(-th) * mp.cos(th)
    print(f"k={k}: M_k / leading law = {mp.nstr(val(tab[k]) / lead, 8)}   (should tend to 1, up to the kappa_1/sqrt(N) term and phase)")
