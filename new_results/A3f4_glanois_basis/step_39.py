import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 1600
cc.N = 1000
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 200
from flint import fmpz, fmpz_mat

def mpc_acb(z): return mp.mpc(mp.mpf(z.real.mid().str(200,radius=False,more=True)), mp.mpf(z.imag.mid().str(200,radius=False,more=True)))

def cmzv_word(xs, etas):
    p = len(xs); run = 0; signs=[0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

# B^2 weight-3 basis: zeta(3; -1), zeta(1,2; 1,-1), zeta(2,1; 1,-1)
w3 = cmzv_word((3,),(4,)); w12 = cmzv_word((1,2),(0,4)); w21 = cmzv_word((2,1),(0,4))
print('words:', w3, w12, w21)
vals = [mpc_acb(cc.h(w)) for w in (w3, w12, w21)]
z3 = mpc_acb(cc.h((-1,-1,0)))
print('zeta(3;-1) value:', mp.nstr(vals[0], 15), ' vs -3/4 zeta3:', mp.nstr(-mp.mpf('0.75')*mp.zeta(3),15))

C = 10**150
M = fmpz_mat(4, 6)
allv = [z3]+vals
for i,v in enumerate(allv):
    M[i,i] = fmpz(1)
    M[i,4] = fmpz(int(mp.nint(C*v.real)))
    M[i,5] = fmpz(int(mp.nint(C*v.imag)))
R = M.lll(delta=0.99)
for r in range(4):
    coeffs = [int(R[r,c]) for c in range(4)]
    res = sum(c*v for c,v in zip(coeffs, allv))
    print(coeffs, 'resid:', mp.nstr(abs(res),3))
