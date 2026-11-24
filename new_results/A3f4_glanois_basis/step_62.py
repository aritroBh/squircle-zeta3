import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 3600
ctx.cap = 5200
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 80
import itertools

def get(w, N):
    cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
    cc.N = N
    return cc.h(w)

def midmag(d):
    re = mp.mpf(d.real.mid().str(12,radius=False)); im = mp.mpf(d.imag.mid().str(12,radius=False))
    return mp.sqrt(re*re+im*im)

def cmzv_word(xs, etas):
    p = len(xs); run = 0; signs=[0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

T2pi = 2*arb.pi()*acb(0,1)
basis = []
for p in (1,2,3):
    for xs in itertools.product([1,2,3], repeat=p):
        if sum(xs)!=3: continue
        for eps in itertools.product([1,-1], repeat=p):
            basis.append((f'z{xs};{eps}', xs, eps, 0))
for p in (1,2):
    for xs in itertools.product([1,2], repeat=p):
        if sum(xs)!=2: continue
        for eps in itertools.product([1,-1], repeat=p):
            basis.append((f'z{xs};{eps}*T', xs, eps, 1))
for eps in itertools.product([1,-1], repeat=1):
    basis.append((f'z(1);{eps}*T^2', (1,), eps, 2))
basis.append(('T^3', (), (), 3))

def basis_word_s(xs, eps):
    etas = [0 if e==1 else 4 for e in eps]
    etas[-1] = (etas[-1]+1)%8
    return cmzv_word(xs, etas)

freewords = [(7,b,c) for b in (4,5,6,7) for c in range(8)]
all_words = [basis_word_s(b[1],b[2]) if b[1] else () for b in basis] + freewords
bad = []
for i,w in enumerate(all_words):
    s = basis[i][3] if i < 27 else 0
    v_hi = get(w, 4200)*T2pi**s
    v_lo = get(w, 2100)*T2pi**s
    dm = midmag(v_hi-v_lo)
    if dm > mp.mpf('1e-300'):
        bad.append((i, w, dm))
print('bad count:', len(bad))
for i,w,dm in bad:
    print(i, w, mp.nstr(dm,3))
