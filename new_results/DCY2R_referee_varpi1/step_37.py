import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
from step_07_partb import momentsum_of, ModelRatio_of, T1ratio_of, D_of, y_of, ell_of, mu_p_of, H_of
PI = arb.pi()
for N in [400, 10000]:
    ms = momentsum_of(N)
    MR = ModelRatio_of(N)
    TR = T1ratio_of(N)
    print(f"N={N}: momentsum = {ms.str(8)}")
    print(f"   ModelRatio = {MR.str(8)}")
    print(f"   T1ratio    = {TR.str(8)}")
    print(f"   D(N)       = {D_of(N).str(8)}")
    # mu_0 .. mu_3
    l2 = ell_of(N,2)
    for p in [0,1,2,3]:
        print(f"   mu_{p} = {mu_p_of(N,p).str(8)}")
