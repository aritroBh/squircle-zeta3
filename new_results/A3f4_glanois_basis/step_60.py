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

def midmag(d):
    re = mp.mpf(d.real.mid().str(15,radius=False)); im = mp.mpf(d.imag.mid().str(15,radius=False))
    return mp.sqrt(re*re+im*im)

for w in [(7,4,0), (-1,-1,5), (5,5,5)]:
    zs = [get(w, N) for N in (1000, 2000, 4000)]
    d1 = midmag(zs[0]-zs[1]); d2 = midmag(zs[1]-zs[2])
    print(w, 'N1000vs2000:', mp.nstr(d1,3), ' N2000vs4000:', mp.nstr(d2,3), ' => digits per 1000 terms ~', mp.nstr(-mp.log10(d2)/4,4) if d2>0 else '>1500')
t=time.time(); get((7,4,0), 4000); print('(7,4,0) N=4000 time:', f'{time.time()-t:.2f}s')
