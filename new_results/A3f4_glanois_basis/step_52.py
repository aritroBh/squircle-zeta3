import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 6000
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 100

vals = {}
for N in (2000, 4000):
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    t=time.time()
    z = cc.h((-1,-1,2))
    vals[N] = z
    print(f'N={N}: time {time.time()-t:.1f}s')
d = vals[2000] - vals[4000]
print('Li3(i) N2000 vs N4000 diff:', abs(d).str(4))
print('expect ~0.855^2000 = 10^%.0f' % (2000*mp.log10(0.855)))
