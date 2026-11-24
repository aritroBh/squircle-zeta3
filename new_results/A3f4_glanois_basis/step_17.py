import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 2400
cc.N = 1200
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 80

for N in (800, 1200, 2000):
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    t = time.time()
    z = cc.h((-1,-1,1))
    err = abs(cc.approx(z) - mp.polylog(3, mp.exp(mp.pi*mp.j/4)))
    print(f'N={N}: radius={z.real.rad().str(3,radius=False)} miderr80d={mp.nstr(err,3)} time={time.time()-t:.1f}s')
