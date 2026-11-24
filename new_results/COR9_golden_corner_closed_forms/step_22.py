import os, time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_16-22: cor9_lib.py written; class_number validated on 62 fields (h(103)=1 confirmed via norm-2 element (477,47)); fund_unit CF-based.\n")
open(WD+'result.md').write('''# COR9 -- number-level closed forms of M3 at the three remaining classical corners + larger census

Work dir `research/08_theory/COR9_golden_corner_closed_forms/`. Labels: **PROVED** = exact algebra;
**NUMERICAL** = lattice relations identical at two precisions with residual at working precision.
Values: HPL = G1's proved B3 term lists by exact series; quad = independent tanh-sinh quadrature (COR7's
precision-boosted integrand). Machinery: cor9_lib.py (reuses cor7_lib/cor7_pt/cor6; adds CF-based
fundamental units, diop_DN norm-l elements, class numbers by form cycles -- validated on 62 fields).

(status: in progress)
''')
print("workspace initialized")

# ---- Task 1, step A: M3, W4, W3, W2 at c0 = (sqrt2-1)^2, dps 1150, via HPL series ----
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 1150
import cor7_lib as L7
t0 = time.time()
rt2 = mp.sqrt(2); b = rt2 - 1; c = b*b   # 3 - 2 sqrt2
print("c0 =", mp.nstr(c, 30))
print("alpha = 2 arccos c0 =", mp.nstr(2*mp.acos(c), 30))
M3, W4, W3, W2 = L7.M3_weights(c)
print(f"HPL series done, {time.time()-t0:.0f}s")
print("M3 =", mp.nstr(M3, 40))
print("W4 =", mp.nstr(W4, 30))
print("W3 =", mp.nstr(W3, 30))
print("W2 =", mp.nstr(W2, 30))
WDp = WD
for nm, v in [("M3", M3), ("W4", W4), ("W3", W3), ("W2", W2)]:
    open(WDp+f"b2sq_{nm}_hpl_1150.txt", "w").write(mp.nstr(v, 1145))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_23: b2sq M3/W4/W3/W2 at dps 1150 via HPL ({time.time()-t0:.0f}s).\n")
print("saved")
