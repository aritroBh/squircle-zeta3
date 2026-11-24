"""w6 identification step 4: c_3, c_2, c_1, c_0 against the pure level-2 weight-6 basis (13 elements, plus variants)
at two precisions (D=200 with the dps-250 run; D=185 with the dps-220 run).  Printed numbers come from this script."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
for j in (3, 2, 1, 0):
    for run, D in ((250, 200), (220, 185)):
        mp.mp.dps = D + 15          # after imports
        C = consts()
        B = basis13(C)
        cj = mp.mpf(coef(run, j))
        mp.mp.dps = D + 15
        print(f"c_{j}: run {run}, D={D}, c = {mp.nstr(cj, 30)}", flush=True)
        rel, res, tol = pslq_try(cj, B, D)
        show("13-basis", B, rel, res, tol, D)
