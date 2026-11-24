
import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import E_series, tompf
flint.ctx.cap = 64
ctx.dps = 110; mp.mp.dps = 110
N = 31
t0=time.time()
E2f, _ = E_series(1, 0, N, h=arb(1)/64, T=4.8)
E3f, _ = E_series(2, 0, N, h=arb(1)/64, T=4.8)
E2g, _ = E_series(1, 0, N, h=arb(1)/32, T=4.8)
E3g, _ = E_series(2, 0, N, h=arb(1)/32, T=4.8)
t1=time.time()
print(f"4 runs in {t1-t0:.0f}s")
E2=[tompf(x,120) for x in E2f]; E3=[tompf(x,120) for x in E3f]
E2b=[tompf(x,120) for x in E2g]; E3b=[tompf(x,120) for x in E3g]
for nm, A, B in [('E2(M_k)', E2, E2b), ('E3(edge)', E3, E3b)]:
    rel = [abs(A[k]-B[k])/max(abs(A[k]),mp.mpf('1e-300')) for k in range(N)]
    kmax = max(range(1,N), key=lambda k: rel[k])
    print(f"{nm}: worst rel h-diff = {mp.nstr(rel[kmax],3)} at k={kmax}; |A_30|={mp.nstr(abs(A[30]),4)}")
    # radius (ball) growth report
json.dump({'E2':[str(x) for x in E2], 'E3':[str(x) for x in E3]}, open('EV_coeffs_K30.json','w'))
print("root test |M_k|^(1/k) (2D, expect ->1):", [mp.nstr(abs(E2[k])**(1/k),5) for k in [6,10,15,20,25,30]])
print("root test |E_k|^(1/k) (3D edge):", [mp.nstr(abs(E3[k])**(1/k),5) for k in [6,10,15,20,25,30]])
