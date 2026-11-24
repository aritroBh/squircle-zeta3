
import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import E_series, tompf

ctx.dps = 110; mp.mp.dps = 110
N = 31
t0 = time.time()
E2f, _ = E_series(1, 0, N, h=arb(1)/64, T=4.8)   # 2D squircle E2 (M_k)
t1 = time.time()
E3f, _ = E_series(2, 0, N, h=arb(1)/64, T=4.8)   # 3D edge E
t2 = time.time()
E2f2, _ = E_series(1, 0, N, h=arb(1)/32, T=4.8)
E3f2, _ = E_series(2, 0, N, h=arb(1)/32, T=4.8)
t3 = time.time()
E2 = [tompf(x, 120) for x in E2f]; E3 = [tompf(x, 120) for x in E3f]
E2b = [tompf(x, 120) for x in E2f2]; E3b = [tompf(x, 120) for x in E3f2]
print(f"times: E2 {t1-t0:.0f}s, E3 {t2-t1:.0f}s, re-run h=1/32 {t3-t2:.0f}s")
# h-convergence
for nm, A, B in [('E2(M_k)', E2, E2b), ('E3(edge)', E3, E3b)]:
    diffs = [abs(A[k]-B[k]) for k in range(N)]
    rel = [diffs[k]/max(abs(A[k]), mp.mpf('1e-90')) for k in range(N)]
    print(f"{nm}: worst rel h-diff = {mp.nstr(max(rel[1:]),3)} at k={rel[1:].index(max(rel[1:]))+1}")
json.dump({'E2':[str(x) for x in E2], 'E3':[str(x) for x in E3]}, open('EV_coeffs_K30.json','w'))
# quick root test on M_k
print("root test |M_k|^(1/k) (2D, should -> 1):")
for k in [5,10,15,20,25,30]:
    print(f"  k={k}: |M_k|={mp.nstr(abs(E2[k]),6)}  root={mp.nstr(abs(E2[k])**(1/k),8)}")
print("root test |E3_k|^(1/k) (3D edge):")
for k in [5,10,15,20,25,30]:
    print(f"  k={k}: |E_k|={mp.nstr(abs(E3[k]),6)}  root={mp.nstr(abs(E3[k])**(1/k),8)}")
