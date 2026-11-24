import sys, os, pickle
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
mp.mp.dps = 300
from flint import fmpz, fmpz_mat
import itertools

def cmzv_word(xs, etas):
    p = len(xs); run = 0; signs=[0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

def mpc_acb(z): return mp.mpc(mp.mpf(z.real.mid().str(300,radius=False,more=True)), mp.mpf(z.imag.mid().str(300,radius=False,more=True)))

T2pi = 2*arb.pi()*acb(0,1)
# B^4 basis at weight 3: zeta(x; 1,...,1, xi4)*(2 pi i)^s, xi4 = i = letter 2
B4 = []
for p in (1,2,3):
    for xs in itertools.product([1,2,3], repeat=p):
        if sum(xs)!=3: continue
        etas = [0]*(p-1)+[2]
        B4.append((f'z{xs};xi4', cmzv_word(xs, etas), 0))
for p in (1,2):
    for xs in itertools.product([1,2], repeat=p):
        if sum(xs)!=2: continue
        etas = [0]*(p-1)+[2]
        B4.append((f'z{xs};xi4*T', cmzv_word(xs, etas), 1))
B4.append(('z(1);xi4*T^2', cmzv_word((1,), (2,)), 2))
B4.append(('T^3', (), 3))
print('B4 size:', len(B4))
B4vals = []
for label, w, s in B4:
    v = cc.h(w) if w else acb(1)
    v = v * T2pi**s
    B4vals.append(mpc_acb(v))
    print(label, w)

C = 10**250
def lll_relations(vals):
    m = len(vals)
    M = fmpz_mat(m, m+2)
    for i,v in enumerate(vals):
        M[i,i] = fmpz(1)
        M[i,m] = fmpz(int(mp.nint(C*v.real)))
        M[i,m+1] = fmpz(int(mp.nint(C*v.imag)))
    R = M.lll(delta=0.99)
    rows = []
    for r in range(m):
        coeffs = [int(R[r,c]) for c in range(m)]
        tail = max(abs(int(R[r,m])), abs(int(R[r,m+1])))
        rows.append((coeffs, tail))
    return rows

z3 = mpc_acb(cc.h((-1,-1,0)))
rows = lll_relations([z3]+B4vals)
rows.sort(key=lambda rt: rt[1])
for coeffs, tail in rows[:3]:
    res = sum(c*v for c,v in zip(coeffs, [z3]+B4vals))
    print('height=', max(abs(c) for c in coeffs), 'tail=', tail, 'resid=', mp.nstr(abs(res),3))
    print('   coeffs:', coeffs)
