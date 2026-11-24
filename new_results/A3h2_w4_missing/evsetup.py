"""Shared evaluator setup for A3h2. Import AFTER setting desired prec/N via setup(prec,cap,N)."""
import sys, os
WORK='<project>/research/08_theory/A3h2_w4_missing'
BASE='<project>'
os.chdir(WORK)
sys.path.insert(0, BASE+'/paper/squircle_note/scripts')
sys.path.insert(0, BASE+'/research/08_theory/A3h_weight4_target/Q')
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp

def setup(prec=700, cap=1000, N=800):
    ctx.prec=prec; ctx.cap=cap; cc.N=N
    cc.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
    cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
    for f in (cc.series,cc.half,cc.g_at_one,cc.h):
        f.cache_clear()

def a2mp(z, digits=100):
    return mp.mpc(mp.mpf(z.real.mid().str(digits,radius=False,more=True)),
                  mp.mpf(z.imag.mid().str(digits,radius=False,more=True)))

GAL={1:1,3:-1,5:-1,7:1}   # sign of sqrt2 under sigma_k
def sig(w,k):
    """Galois sigma_k on a word: letter j -> k*j mod 8, -1 fixed."""
    return tuple((-1 if x==-1 else (k*x)%8) for x in w)

def h(w):
    return cc.h(tuple(w))

def h_k(w,k):
    """evaluate h at Galois embedding k = h(sigma_k(w)) with standard roots."""
    return cc.h(sig(w,k))

def log(step,msg):
    s=f'{step}: {msg}'
    print(s,flush=True)
    with open(WORK+'/progress.log','a') as f: f.write(s+'\n')
