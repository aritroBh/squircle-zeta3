import sys, os, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 1200
import sympy as sp
import itertools

def mpc_acb(z): return mp.mpc(mp.mpf(z.real.mid().str(1100,radius=False,more=True)), mp.mpf(z.imag.mid().str(1100,radius=False,more=True)))
def cmzv_word(xs, etas):
    p = len(xs); run = 0; signs=[0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

T2pi = 2*arb.pi()*acb(0,1)
B4 = []
for p in (1,2,3):
    for xs in itertools.product([1,2,3], repeat=p):
        if sum(xs)!=3: continue
        B4.append((f'z{xs};xi4', cmzv_word(xs, [0]*(p-1)+[2]), 0))
for p in (1,2):
    for xs in itertools.product([1,2], repeat=p):
        if sum(xs)!=2: continue
        B4.append((f'z{xs};xi4*T', cmzv_word(xs, [0]*(p-1)+[2]), 1))
B4.append(('z(1);xi4*T^2', cmzv_word((1,), (2,)), 2))
B4.append(('T^3', (), 3))
B4vals = []
for label, w, s in B4:
    v = cc.h(w) if w else acb(1)
    B4vals.append(mpc_acb(v*T2pi**s))
z3 = mpc_acb(cc.h((-1,-1,0)))

C = 10**500
m = 9
rows = []
for i,v in enumerate([z3]+B4vals):
    r = [0]*(m+2)
    r[i] = 1
    r[m] = int(mp.nint(C*v.real))
    r[m+1] = int(mp.nint(C*v.imag))
    rows.append(r)
t=time.time()
R = sp.Matrix(rows).lll(delta=sp.Rational(99,100))
print('lll time', time.time()-t)
for r in R.tolist()[:3]:
    res = sum(c*v for c,v in zip(r, [z3]+B4vals))
    h = max(abs(c) for c in r)
    print('height 10^%.1f' % mp.log10(h), 'residual:', mp.nstr(abs(res),3), 'ratio 10^%.1f' % mp.log10(abs(res)/h))
