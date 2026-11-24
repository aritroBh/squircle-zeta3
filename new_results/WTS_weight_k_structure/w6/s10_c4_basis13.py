"""w6 identification step 3: c_4 against the full 13-element weight-6 alternating-MZV basis, at two precisions
(D=200 with c_4 from the independent dps-250 run; D=185 with c_4 from the independent dps-220 run).
Also a check at D=380 with the dps-400 run on the 8-element subset.  All printed numbers come from this script."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
for run, D in ((250, 200), (220, 185)):
    mp.mp.dps = D + 15          # after imports
    C = consts()
    B = basis13(C)
    c4 = mp.mpf(coef(run, 4))
    mp.mp.dps = D + 15
    print(f"run {run}: c_4 =", mp.nstr(c4, 40), " mp.mp.dps =", mp.mp.dps)
    rel, res, tol = pslq_try(c4, B, D)
    show("13-basis", B, rel, res, tol, D)
# closed form residual at several precisions
for run, D in ((400, 380), (250, 235), (220, 205)):
    mp.mp.dps = D
    lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2
    cf = (168 * lam ** 6 + 420 * pi2 * lam ** 4 + 378 * pi2 ** 2 * lam ** 2 + 61 * pi2 ** 3 + 5040 * z3 * lam ** 3
          + 2520 * pi2 * z3 * lam + 3780 * z3 ** 2 + 22680 * z5 * lam) / 362880
    print(f"closed form - c_4(run {run}) = {mp.nstr(cf - mp.mpf(coef(run, 4)), 5)}   (mp.mp.dps = {mp.mp.dps})")
