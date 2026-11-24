import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 5000
ctx.cap = 7000
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 80

def get(w, N):
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    return cc.h(w)

for w in [(7,4,0), (-1,-1,5), (5,5,5)]:
    zs = [get(w, N) for N in (1000, 2000, 4000)]
    d1 = abs(zs[0]-zs[1]); d2 = abs(zs[1]-zs[2])
    r = mp.log10(mp.mpf(d1.str(5))) / mp.log10(mp.mpf(d2.str(5))) if d2>0 else mp.inf
    print(w, 'N1000vs2000:', d1.str(3), ' N2000vs4000:', d2.str(3), ' digits/1000 terms ~', mp.nstr(-mp.log10(mp.mpf(d2.str(5)))/4,5) if d2>0 else 'inf')
print('time check:')
t=time.time(); get((7,4,0), 4000); print('(7,4,0) N=4000:', f'{time.time()-t:.1f}s')
