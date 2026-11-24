import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import V_series, E_series, tompf
flint.ctx.cap = 64
ctx.dps = 110; mp.mp.dps = 110
os.chdir('<project>/research/08_theory/SRF2_cube_residue')
N = 41
t0 = time.time()
Vf, Vh = V_series(0, N, h=arb(1)/64, T=4.6)
t1 = time.time()
Vg, _ = V_series(0, N, h=arb(1)/32, T=4.6)
t2 = time.time()
E3f, _ = E_series(2, 0, N, h=arb(1)/64, T=4.6)
t3 = time.time()
print(f"V h=1/64: {t1-t0:.0f}s, V h=1/32: {t2-t1:.0f}s, E3: {t3-t2:.0f}s")
V = [tompf(x,120) for x in Vf]; Vb = [tompf(x,120) for x in Vg]
E3 = [tompf(x,120) for x in E3f]
rel = [abs(V[k]-Vb[k])/max(abs(V[k]),mp.mpf('1e-300')) for k in range(N)]
print("V worst rel h-diff per decade:")
for k in [10, 20, 30, 40]:
    print(f"  k={k}: {mp.nstr(rel[k],3)}")
json.dump({'V':[str(x) for x in V], 'E3':[str(x) for x in E3]}, open('EV_coeffs_K40.json','w'))
print("saved EV_coeffs_K40.json")
