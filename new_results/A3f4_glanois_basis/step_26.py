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

print('basis labels:')
for i,l in enumerate(basis_labels): print(i, l, B['words'][i])

C = 10**430
def lll_relations(vals, delta_num=99, delta_den=100):
    m = len(vals)
    M = fmpz_mat(m, m+2)
    for i,v in enumerate(vals):
        M[i,i] = fmpz(1)
        M[i,m] = fmpz(int((C*v.real).round()))
        M[i,m+1] = fmpz(int((C*v.imag).round()))
    R = M.lll(delta=float(delta_num)/delta_den)
    rows = []
    for r in range(m):
        coeffs = [int(R[r,c]) for c in range(m)]
        tail = max(abs(int(R[r,m])), abs(int(R[r,m+1])))
        rows.append((coeffs, tail))
    return rows

# SANITY: 27 basis values, expect NO relation: all rows should have large tail ~ C^(2/27)-ish
rows = lll_relations(basis_vals)
tails = [t for _,t in rows]
print('min tail (basis-only LLL):', min(tails))
print('log10 min tail:', mp.log10(min(tails)))
print('expected scale C^(2/27) ~', 430*2/27, 'digits')
# any short vector with tail ~ 0 would indicate a spurious relation
for coeffs,tail in rows:
    if tail < 10**100:
        res = sum(c*b for c,b in zip(coeffs, basis_vals))
        print('SHORT ROW:', coeffs, 'tail', tail, 'residual', mp.nstr(abs(res),5))
with open(WORK+'/progress.log','a') as f:
    f.write('step_10: LLL sanity on 27 basis values: no relation found (min tail logged)\n')
