import sys, os, time
WORK='<project>/research/08_theory/A3h2_w4_missing'
os.chdir(WORK)
BASE='<project>'
sys.path.insert(0, BASE+'/paper/squircle_note/scripts')
QH=BASE+'/research/08_theory/A3h_weight4_target/Q'
sys.path.insert(0, QH)
t0=time.time()
import check_conversion as cc
from flint import ctx, acb, arb
print('import ok', time.time()-t0)
# reset traps
ctx.prec=700; ctx.cap=1000; cc.N=800
cc.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
for f in (cc.series,cc.half,cc.g_at_one,cc.h):
    f.cache_clear()
import mpmath as mp
mp.mp.dps=60
t0=time.time()
z=cc.h((7,4,0))
print('h(7,4,0) at N=800:', time.time()-t0,'s')
print(mp.nstr(cc.approx(z),30))
t0=time.time()
z2=cc.h((7,4,7,7))
print('h(7,4,7,7):', time.time()-t0,'s')
print(mp.nstr(cc.approx(z2),30))
# independent check: zeta(3)=h(-1,-1,0)
t0=time.time()
z3=cc.h((-1,-1,0))
print('h(-1,-1,0) vs zeta(3):', mp.nstr(abs(cc.approx(z3)-mp.zeta(3)),5), time.time()-t0)
