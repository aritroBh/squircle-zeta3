import os, pickle
from fractions import Fraction as Fr
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][:27]]
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']
with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)
with open(WORK+'/pivots.pickle','rb') as f: P = pickle.load(f)
rows_int = P['rows_int']; pivot = P['pivot_rows']; nonp = P['nonp']

# (1) verify expansion of (7,4,1)
w = (7,4,1)
i = freewords.index(w)
res = free_vals[i] - sum(q*b for q,b in zip(E[w]['coeffs'], basis_vals))
print('expansion (7,4,1) residual:', mp.nstr(abs(res),3))
# and the integer LLL relation form: a0*w + sum a_j b_j = 0
a0 = E[w]['a0']
res2 = a0*free_vals[i] + sum(int(-q*a0)*b for q,b in zip(E[w]['coeffs'], basis_vals))
print('integer form residual:', mp.nstr(abs(res2),3))
# rows_int should equal [-q*a0]
print('rows_int matches:', all(int(-q*a0)==x for q,x in zip(E[w]['coeffs'], rows_int[i])))

# (2) verify row relation arithmetic exactly
def solve27(A, b):
    n = 27
    M = [row[:] + [b[i]] for i,row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        inv = 1/M[col][col]
        M[col] = [x*inv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f*bb for a,bb in zip(M[r], M[col])]
    return [M[i][n] for i in range(n)]

Pmat = [[Fr(x) for x in rows_int[j]] for j in pivot]
k = nonp[0]
c = solve27([row[:] for row in Pmat], [Fr(x) for x in rows_int[k]])
diff = [rows_int[k][j] - sum(ci*rows_int[pi][j] for ci,pi in zip(c,pivot)) for j in range(27)]
print('row relation arithmetic exact?', all(d==0 for d in diff))
