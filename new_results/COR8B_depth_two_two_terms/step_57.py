
import sys, pickle
from fractions import Fraction as Fr
br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
# rows = 78 reduced rows (dict). compute left null: f in Q^81 with f.r = 0 for all r
# build dense matrix of pivot rows
piv_rows = {}
for name, row in br:
    piv_rows[min(row)] = (row, name)
piv_cols = sorted(piv_rows)
print("pivot cols:", piv_cols)
# a functional f: values on pivot cols determined by values on free cols
free_cols = [j for j in range(81) if j not in piv_rows]
print("free cols:", free_cols)
# for each free col j, functional f_j: f_j[j]=1, f_j[piv] = -row[piv coeff...]
# constraint: sum_i f_i * row_i = 0 for each row: f[min] * row[min] + sum_{other} f[other]*row[other] = 0
nulls = []
for j in free_cols:
    f = {j: Fr(1)}
    for pc in piv_cols:
        row, name = piv_rows[pc]
        # f·row = 0 -> f[pc]*row[pc] = -sum_{j' != pc} f[j']*row[j']
        s = sum(f.get(jj,Fr(0))*x for jj,x in row.items() if jj != pc)
        f[pc] = -s/row[pc]
    nulls.append(f)
for idx,f in enumerate(nulls):
    nz = {i:x for i,x in f.items() if x}
    print(f"null {idx}: ", {i:str(x) for i,x in sorted(nz.items())})
# apply to targets
T = pickle.load(open('targets81.pkl','rb'))
for e in (-1,0,1):
    vals = [sum(f.get(i,Fr(0))*x for i,x in enumerate(T[e])) for f in nulls]
    print("e=",e,"null evaluations:", [str(v) for v in vals])
