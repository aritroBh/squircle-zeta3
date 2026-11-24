import sys, os, pickle
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
mp.mp.dps = 620
from flint import fmpz, fmpz_mat

with open(WORK+'/basis_vals.pickle','rb') as f: B = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in B['vals']]

def mpc_acb(z): return mp.mpc(mp.mpf(z.real.mid().str(600,radius=False,more=True)), mp.mpf(z.imag.mid().str(600,radius=False,more=True)))

targets = {
  'zeta3=(-1,-1,0)': mpc_acb(cc.h((-1,-1,0))),
  'Li3(i)=(-1,-1,2)': mpc_acb(cc.h((-1,-1,2))),
  'Li3(xi8^3)=(-1,-1,3)': mpc_acb(cc.h((-1,-1,3))),
}

C = 10**430
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

for name, tv in targets.items():
    vals = [tv] + basis_vals
    rows = lll_relations(vals)
    rows.sort(key=lambda rt: rt[1])
    coeffs, tail = rows[0]
    res = sum(c*v for c,v in zip(coeffs, vals))
    print(name, 'c0=', coeffs[0], 'height=', max(abs(c) for c in coeffs), 'resid=', mp.nstr(abs(res),3))
