import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp
mp.mp.dps = 120

def fresh_roots():
    cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
    cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)

def get(w, prec, N):
    ctx.prec = prec; fresh_roots()
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    return cc.h(w)

for prec, N in ((4200,1400),(4200,3000),(8000,1000),(8000,3000),(8000,6000)):
    li = get((-1,-1,2), prec, N)
    im = li.imag.mid()
    pim = arb.pi()**3/32
    ctx.prec = 800
    diff = abs(im - pim)
    mp.mp.dps=120
    print(f'prec={prec} N={N}: |ImLi3i - pi^3/32| = {mp.nstr(mp.mpf(diff.str(30,radius=False)),4)}')
