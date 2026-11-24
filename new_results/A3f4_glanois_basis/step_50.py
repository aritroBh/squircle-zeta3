import sys, os, itertools, time
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
mp.mp.dps = 300
import sympy as sp

def hmp(w, digs=280):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(digs,radius=False,more=True)), mp.mpf(hv.imag.mid().str(digs,radius=False,more=True)))

T2pi = 2*arb.pi()*acb(0,1)
def T2pi_mp(digs=280):
    return mp.mpc(0, 2*mp.pi)

elems = {
 'Li3(i)': hmp((-1,-1,2)),
 'Li12(1,i)': hmp((-1,2,2)),
 'Li21(1,i)': hmp((2,-1,2)),
 'Li111(1,1,i)': hmp((2,2,2)),
 'Li2(i)T': hmp((-1,2)) * T2pi_mp(),
 'Li11(1,i)T': hmp((2,2)) * T2pi_mp(),
 'Li1(i)T2': hmp((2,)) * T2pi_mp()**2,
 'T3': T2pi_mp()**3,
}
z3 = hmp((-1,-1,0))
C = 10**250
def lll_rel(vals):
    m = len(vals)
    rows = []
    for i,v in enumerate(vals):
        r = [0]*(m+2); r[i]=1
        r[m] = int(mp.nint(C*v.real)); r[m+1] = int(mp.nint(C*v.imag))
        rows.append(r)
    R = sp.Matrix(rows).lll(delta=sp.Rational(99,100))
    out = []
    for r in R.tolist():
        res = sum(c*v for c,v in zip(r, vals))
        hgt = max(abs(c) for c in r)
        out.append((r, hgt, abs(res)))
    return out

names = list(elems)
# try subsets of size 2..6 (of basis elements) + zeta3
best = []
for k in range(2, 7):
    for sub in itertools.combinations(range(8), k):
        vals = [z3] + [elems[names[i]] for i in sub]
        for r, hgt, res in lll_rel(vals):
            if hgt > 0 and res/hgt < mp.mpf('1e-230'):
                best.append((hgt, sub, r, res))
    if best:
        best.sort()
        print('k=',k,'found', len(best))
        break
if best:
    hgt, sub, r, res = best[0]
    print('subset:', [names[i] for i in sub])
    print('coeffs [zeta3]+sub:', r)
    print('height:', hgt, 'residual:', mp.nstr(res,3))
else:
    print('no relation found in any subset up to size 6')
