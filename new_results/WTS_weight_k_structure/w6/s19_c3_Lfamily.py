"""w6 identification step 12: c_3 against the L-family { L^m * (level-2 weight 6-m basis) }, L = log(1+sqrt2),
33 elements (no sqrt2) and 66 elements (with sqrt2 multiples), by LLL (python-flint) at two precisions:
(run 400, D=250) and (run 250, D=225)  [S6 cache has 275 digits].  A relation is reported with its max coefficient; with n elements and D digits
a spurious relation has max coefficient ~10^(D/n).  Printed numbers from this script."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel
for sq in (False, True):
    for run, D in ((400, 250), (250, 225)):
        mp.mp.dps = D + 20
        C = consts()
        Lval = mp.log(1 + mp.sqrt(2))
        B = Lfamily(C, Lval, sqrt2=sq)
        c3 = mp.mpf(coef(run, 3))
        mp.mp.dps = D + 20
        n = len(B) + 1
        print(f"c_3, run {run}, D={D}, {n} numbers (sqrt2 multiples: {sq}), mp.mp.dps = {mp.mp.dps}; spurious-relation size ~10^{D/n:.1f}", flush=True)
        R = find_rel(c3, B, D)
        if not R:
            print("   no relation with residual < 1e-%d" % (D - 25))
        for r in R[:3]:
            print(f"   relation: max|coeff| = {r['maxcoeff']}, residual {mp.nstr(r['res'],3)} (dps {mp.mp.dps}); target coeff {r['rel'][0]}")
            if r["maxcoeff"] < 10 ** 6: print("     ", r["terms"])
