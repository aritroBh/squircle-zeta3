import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
from flint import ctx
ctx.prec = 2000
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
import mpmath as mp
mp.mp.dps = 50

for N in (500, 1000, 1500, 2000):
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    t = time.time()
    z = cc.h((-1,-1,1))
    err = abs(cc.approx(z) - mp.polylog(3, mp.exp(mp.pi*mp.j/4)))
    print(f'N={N}: radius_re={z.real.rad().str(3,radius=False)} miderr={mp.nstr(err,3)} time={time.time()-t:.1f}s')
