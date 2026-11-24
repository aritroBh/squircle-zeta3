# Re-run (not independent) of the note's exact reducer at k = 4, 5; read-only import, no bytecode written.
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from ibp_reduce import reduce_k
for k in [6]:
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
