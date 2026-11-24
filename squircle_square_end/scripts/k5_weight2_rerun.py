"""k5_weight2_rerun.py -- re-run (not independent) of the exact reducer ibp_reduce.reduce_k at k = 4, 5.

Paper: "Investigating zeta(3)", Proposition 1 (weight two, k <= 5), Table III entry "exact re-run (k=4,5)".
For each k, the weight-2 piece of the reduced M_k integrand is compared with the kernel of
T = (pi^2 + 6 log^2 2)/24; the difference must simplify to exactly 0.
Run from this folder:  python k5_weight2_rerun.py      (about 2 minutes; reference output in
../outputs/k5_weight2_rerun.out)
"""
import os
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from ibp_reduce import reduce_k

for k in [4, 5]:
    res = reduce_k(k, verbose=False)
    print("k =", k, "keys:", sorted(res.keys()), flush=True)
    for key, val in res.items():
        if sum(key) != 1:
            continue
        (v,) = val.free_symbols
        S2 = sp.sqrt(2); a = 1 + S2; b = S2 - 1
        Rw2 = 1/(2*v) + 1/(2*(v - a)) - 1/(v + a) - 1/(2*(v + b))
        d = sp.simplify(sp.radsimp(sp.together(val - Rw2)))
        print("  weight-2 key", key, " minus T-kernel ->", d, "  numeric at v=0.7, 1.9:",
              [sp.N((val - Rw2).subs(v, x), 30) for x in (sp.Rational(7, 10), sp.Rational(19, 10))], flush=True)
    print("  degree-0 (weight-1) pieces:", {kk: vv for kk, vv in res.items() if sum(kk) == 0}, flush=True)
