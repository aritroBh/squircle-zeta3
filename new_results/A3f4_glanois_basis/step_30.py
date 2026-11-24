import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 620
from flint import fmpz, fmpz_mat
from fractions import Fraction as Fr

with open(WORK+'/basis_vals.pickle','rb') as f: B = pickle.load(f)
with open(WORK+'/free_vals_hi.pickle','rb') as f: F = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in B['vals']]
free_vals = {w: mp.mpc(mp.mpf(re), mp.mpf(im)) for w,(re,im) in F['freevals'].items()}

C = 10**550
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

w = (7,5,2)
vals = [free_vals[w]] + basis_vals
rows = lll_relations(vals)
rows.sort(key=lambda rt: rt[1])
for coeffs, tail in rows[:6]:
    print('c0=', coeffs[0], 'height=', max(abs(c) for c in coeffs), 'tail log10=', mp.nstr(mp.log10(tail),5) if tail else '-inf')
    res = sum(c*v for c,v in zip(coeffs, vals))
    print('   residual=', mp.nstr(abs(res),4))
with open(WORK+'/progress.log','a') as f:
    f.write('step_14: retry (7,5,2) with C=1e550\n')
