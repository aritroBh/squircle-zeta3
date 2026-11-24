
import sys, pickle
from fractions import Fraction as Fr
br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
piv_rows = {}
for name, row in br:
    piv_rows[min(row)] = (row, name)
piv_cols = sorted(piv_rows)
free_cols = [j for j in range(81) if j not in piv_rows]
nulls = []
for j in free_cols:
    f = {j: Fr(1)}
    for pc in reversed(piv_cols):      # decreasing order!
        row, name = piv_rows[pc]
        s = sum(f.get(jj,Fr(0))*x for jj,x in row.items() if jj != pc)
        f[pc] = -s/row[pc]
    nulls.append(f)
# sanity: check f·row = 0 for all rows
bad = 0
for f in nulls:
    for name, row in br:
        s = sum(f.get(j,Fr(0))*x for j,x in row.items())
        if s != 0: bad += 1
print("null functional check failures:", bad, "(should be 0)")
T = pickle.load(open('targets81.pkl','rb'))
for e in (-1,0,1):
    vals = [sum(f.get(i,Fr(0))*x for i,x in enumerate(T[e])) for f in nulls]
    print("e=",e,"null evaluations:", [str(v) for v in vals])
