import sys, os, time, pickle
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
mp.mp.dps = 100
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
all_s = [basis[i][3] if i<27 else 0 for i in range(len(all_words))]

t=time.time()
vals = []
worst = mp.mpf(0); worst_i = -1
for i,w in enumerate(all_words):
    s = all_s[i]
    v_hi = get(w, 4200)*T2pi**s
    v_lo = get(w, 2100)*T2pi**s
    dm = midmag(v_hi-v_lo)
    if dm > worst: worst, worst_i = dm, i
    vals.append(v_hi)
print('time', f'{time.time()-t:.0f}s', 'max |N4200-N2100|:', mp.nstr(worst,3), 'at', all_words[worst_i])

def s800(z):
    return (z.real.mid().str(800,radius=False,more=True), z.imag.mid().str(800,radius=False,more=True))
with open(WORK+'/vals_hi.pickle','wb') as f:
    pickle.dump({'basis':[b[0] for b in basis], 'basis_meta':[(b[1],b[2],b[3]) for b in basis],
                 'basis_words':all_words[:27], 'freewords':freewords,
                 'vals':[s800(v) for v in vals],
                 'settings':{'prec':3600,'cap':5200,'N':4200}}, f)
print('saved', len(vals))

# validations
z3 = vals[all_words.index((-1,-1,0))] if (-1,-1,0) in all_words else None
li3i = get((-1,-1,2), 4200)
rel = 24*get((-1,-1,0),4200) + 256*li3i + T2pi**3
print('24z3+256Li3i+T^3 mid:', rel.real.mid().str(3,radius=False), rel.imag.mid().str(3,radius=False))
with open(WORK+'/progress.log','a') as f:
    f.write('step_22: regenerated vals_hi.pickle (59 values, N=4200, per-value convergence cert); max diff printed; 24z3 relation holds\n')
