
from flint import ctx, arb, acb
ctx.dps = 40
from step_07_partb import *
for p in [0, 1, 2, 3, 4, 6, 9]:
    v = mu_p_of(400, p)
    print(f"mu({p}) = {v.str(8)}")
ms = momentsum_of(400)
print("momentsum =", ms.str(8))
print("ModelRatio =", ModelRatio_of(400).str(10))
print("T1ratio =", T1ratio_of(400).str(10))
