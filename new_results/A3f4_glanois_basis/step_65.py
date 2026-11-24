import sys, os, pickle, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp
from fractions import Fraction as Fr

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][:27]]
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']

C = 10**700
def lll_rows(vals):
    m = len(vals)
    rows = []
    for i,v in enumerate(vals):
        r = [0]*(m+2); r[i]=1
        r[m] = int(mp.nint(C*v.real)); r[m+1] = int(mp.nint(C*v.imag))
        rows.append(r)
    return sp.Matrix(rows).lll(delta=sp.Rational(99,100)).tolist()

expansions = {}
t0=time.time()
for w, wv in zip(freewords, free_vals):
    vals = [wv] + basis_vals
    R = lll_rows(vals)
    best = None
    for r in R:
        if r[0] == 0: continue
        tail = max(abs(r[28]), abs(r[29]))
        if best is None or tail < best[1]:
            best = (r, tail)
    r, tail = best
    res = sum(c*v for c,v in zip(r[:28], vals))
    a0 = r[0]
    expr = [Fr(-c, a0) for c in r[1:28]]
    expansions[w] = {'coeffs': expr, 'a0': a0, 'residual': abs(res), 'tail': tail, 'height': max(abs(c) for c in r[:28])}
    print(w, 'a0=', a0, 'H=', max(abs(c) for c in r[:28]), 'tail=10^%.1f'%mp.log10(tail) if tail else 0, 'resid=', mp.nstr(abs(res),3), f'{time.time()-t0:.0f}s')
with open(WORK+'/expansions.pickle','wb') as f:
    pickle.dump(expansions, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_23: LLL expansions of 32 free words in B^8 (sympy exact LLL, C=1e700, values 800 digits)\n')
