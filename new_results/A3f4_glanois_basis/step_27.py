import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 470
from flint import fmpz, fmpz_mat

with open(WORK+'/basis_vals.pickle','rb') as f: B = pickle.load(f)
with open(WORK+'/free_vals_hi.pickle','rb') as f: F = pickle.load(f)
basis_labels = B['basis']; basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in B['vals']]
freewords = F['freewords']; free_vals = {w: mp.mpc(mp.mpf(re), mp.mpf(im)) for w,(re,im) in F['freevals'].items()}

C = fmpz(10)**430
def lll_relations(vals):
    m = len(vals)
    M = fmpz_mat(m, m+2)
    for i,v in enumerate(vals):
        M[i,i] = fmpz(1)
        M[i,m] = fmpz(mp.nint(C*v.real))
        M[i,m+1] = fmpz(mp.nint(C*v.imag))
    R = M.lll(delta=0.99)
    rows = []
    for r in range(m):
        coeffs = [int(R[r,c]) for c in range(m)]
        tail = max(abs(int(R[r,m])), abs(int(R[r,m+1])))
        rows.append((coeffs, tail))
    return rows

rows = lll_relations(basis_vals)
tails = sorted(t for _,t in rows)
print('smallest tails (basis-only LLL):')
for t in tails[:5]: print(' ', t, 'log10=', mp.nstr(mp.log10(t),6) if t>0 else '-inf')
print('largest:', tails[-1])
print('C^(2/27) has', 430*2/27, 'digits; any true relation of height H gives tail ~ H*sqrt(27)')
with open(WORK+'/progress.log','a') as f:
    f.write('step_11: basis-only LLL done; min tail printed (sanity: no relation to height ~1e28)\n')
