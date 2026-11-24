import sys, os, pickle
from math import gcd
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)

with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)
freewords = [(7,b,c) for b in (4,5,6,7) for c in range(8)]

rows_int = []
for w in freewords:
    a0 = E[w]['a0']
    rows_int.append([int(-q*a0) for q in E[w]['coeffs']])  # integers, = -c_j

p = 2**61 - 1
M = [[x % p for x in row] for row in rows_int]  # 32 x 27
# Gaussian elimination to find independent rows (pivot set)
work = []
pivot_rows = []
col_of = {}
used_cols = set()
for i in range(32):
    r = M[i][:]
    for j in range(27):
        if j in used_cols: continue
        if r[j] % p != 0:
            # eliminate r against existing basis rows
            for (bj, brow) in work:
                f = r[bj] * pow(brow[bj], p-2, p) % p
                r = [(a - f*b) % p for a,b in zip(r, brow)]
            if r[j] % p != 0:
                used_cols.add(j)
                work.append((j, r))
                pivot_rows.append(i)
                break
print('num independent rows:', len(pivot_rows))
print('pivot row indices:', pivot_rows)
# the 5 non-pivot rows:
nonp = [i for i in range(32) if i not in pivot_rows]
print('non-pivot (relation) rows:', nonp, [freewords[i] for i in nonp])
import pickle as pk
with open(WORK+'/pivots.pickle','wb') as f:
    pk.dump({'pivot_rows': pivot_rows, 'nonp': nonp, 'rows_int': rows_int}, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_25: mod-p pivot set of 27 rows found; 5 relation rows identified\n')
