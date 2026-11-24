
import sys, time, json
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import V_series, tompf
flint.ctx.cap = 64
ctx.dps = 110; mp.mp.dps = 110
N = 31
t0=time.time()
Vf, Vh = V_series(0, N, h=arb(1)/64, T=4.6)
t1=time.time()
Vg, _ = V_series(0, N, h=arb(1)/32, T=4.6)
t2=time.time()
print(f"V h=1/64: {t1-t0:.0f}s, h=1/32: {t2-t1:.0f}s")
V=[tompf(x,120) for x in Vf]; Vb=[tompf(x,120) for x in Vg]
rel = [abs(V[k]-Vb[k])/max(abs(V[k]),mp.mpf('1e-300')) for k in range(N)]
kmax = max(range(1,N), key=lambda k: rel[k])
print(f"V: worst rel h-diff = {mp.nstr(rel[kmax],3)} at k={kmax}")
old = json.load(open('EV_coeffs_K30.json'))
old['V'] = [str(x) for x in V]
json.dump(old, open('EV_coeffs_K30.json','w'))
# Q_k = -V-related? Q(eps) = -H3/Gamma(1+2eps); V here is 48 eps^2 V coefficient of S: s_k = 48(E_{k-1}+V_{k-2})
# print V_k and root test
print("V_k:", [mp.nstr(V[k],5) for k in range(0,8)])
print("root test |V_k|^(1/k):", [mp.nstr(abs(V[k])**(1/k),5) for k in [6,10,15,20,25,30]])
# s_k
E2=[mp.mpf(x) for x in old['E2']]; E3=[mp.mpf(x) for x in old['E3']]
s = [mp.mpf(24)] + [48*(E3[k-1] + (V[k-2] if k>=2 else 0)) for k in range(1,N)]
print("s_k root test |s_k|^(1/k):", [mp.nstr(abs(s[k])**(1/k),5) for k in [6,10,15,20,25,30]])
claim = {1:'-26.558558386584754433798785874062052986842743535978204080483978',
         2:'-4.6892155184605022499600569810270175081392278040492906513172203885',
         3:'48.287432486014447759215859342074924467116625727590718495681774536',
         4:'-78.06628953961050950490562401928703057980'}
for k in claim: print(f"s{k} diff vs paper: {mp.nstr(mp.mpf(claim[k])-s[k],3)}")
print("s5..s8:", [mp.nstr(s[k],12) for k in range(5,9)])
