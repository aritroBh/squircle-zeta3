import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 300
consts = {k: mp.mpf(v) for k,v in json.load(open('purew3_300.json')).items()}
names = list(consts.keys()); xs = [consts[k] for k in names]
n = len(xs)
from flint import fmpz_mat
B = mp.power(10, 285)
rows=[]
for i in range(n):
    row=[0]*(n+1)
    row[i]=1
    row[-1]=int(mp.nint(B*xs[i]))
    rows.append(row)
M=fmpz_mat(rows)
R = M.lll()
print("nrows,ncols:", R.nrows(), R.ncols())
data=[]
for i in range(R.nrows()):
    v=[R[i,j] for j in range(R.ncols())]
    data.append(v)
data.sort(key=lambda v: abs(v[-1]))
print("smallest |residual*B| rows:")
for v in data[:14]:
    coef=v[:-1]; maxc=max(abs(c) for c in coef)
    print(f"   residB={abs(v[-1]):>28} maxc={maxc:<14} nnz={sum(1 for c in coef if c)}")
log("step: LLL access fixed")
