import sys, os, pickle, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 470
from flint import fmpz, fmpz_mat
from fractions import Fraction as Fr

with open(WORK+'/basis_vals.pickle','rb') as f: B = pickle.load(f)
with open(WORK+'/free_vals_hi.pickle','rb') as f: F = pickle.load(f)
basis_labels = B['basis']; basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in B['vals']]
freewords = F['freewords']; free_vals = {w: mp.mpc(mp.mpf(re), mp.mpf(im)) for w,(re,im) in F['freevals'].items()}

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

expansions = {}   # word -> list of Fractions (coeffs in basis), plus diag info
t0 = time.time()
for w in freewords:
    vals = [free_vals[w]] + basis_vals
    rows = lll_relations(vals)
    # find a row with nonzero w-coeff and tiny tail
    best = None
    for coeffs, tail in rows:
        if coeffs[0] != 0 and tail < 10**30:
            res = sum(c*v for c,v in zip(coeffs, vals))
            key = (abs(coeffs[0]), max(abs(c) for c in coeffs))
            if best is None or key < best[0]:
                best = (key, coeffs, tail, abs(res))
    if best is None:
        print(w, 'NO RELATION FOUND'); continue
    key, coeffs, tail, resabs = best
    a0 = coeffs[0]
    expr = [Fr(-c, a0) for c in coeffs[1:]]
    expansions[w] = {'coeffs': expr, 'a0': a0, 'maxh': max(abs(c) for c in coeffs), 'residual': resabs, 'tail': tail}
    print(w, 'a0=', a0, 'height=', max(abs(c) for c in coeffs), 'residual=', mp.nstr(resabs,4), f'({time.time()-t0:.0f}s)')
with open(WORK+'/expansions.pickle','wb') as f:
    pickle.dump(expansions, f)
print('done', len(expansions))
with open(WORK+'/progress.log','a') as f:
    f.write('step_13: LLL expansions of 32 free words in B^8 basis computed (residuals printed)\n')
