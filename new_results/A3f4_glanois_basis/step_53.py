import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp
mp.mp.dps = 50

def fresh_roots():
    cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
    cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)

def get(w, prec, N):
    ctx.prec = prec
    fresh_roots()
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    return cc.h(w)

ref = get((-1,-1,2), 8000, 6000)
print('error vs N at prec=8000:')
for N in (1000, 1500, 2000, 3000, 4500):
    z = get((-1,-1,2), 8000, N)
    print(f'  N={N}: err {abs(z-ref).str(4)}')
print('error vs prec at N=6000:')
for prec in (3000, 4200, 6000):
    z = get((-1,-1,2), prec, 6000)
    print(f'  prec={prec}: err {abs(z-ref).str(4)}')
