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
    row=[0]*(n+1); row[i]=1; row[-1]=int(mp.nint(B*xs[i])); rows.append(row)
R = fmpz_mat(rows).lll()
data=[]
for i in range(R.nrows()):
    v=[int(R[i,j]) for j in range(R.ncols())]
    data.append(v)
data.sort(key=lambda v: abs(v[-1]))
print("smallest |residual*B| rows (residB, maxc, nnz):")
for v in data[:16]:
    coef=v[:-1]; maxc=max(abs(c) for c in coef); nnz=sum(1 for c in coef if c)
    print("   residB=%.3e maxc=%d nnz=%d"%(abs(v[-1]), maxc, nnz))
# The relation rows have residB ~ 0 (small). Count them
relrows=[v for v in data if abs(v[-1]) < 1e6]
print("\nrows with residB < 1e6 (candidate exact relations):", len(relrows))
log("step: LLL ran; %d candidate relations among 33 pure-w3 consts"%len(relrows))
