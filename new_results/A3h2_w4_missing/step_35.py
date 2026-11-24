import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=180
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
fw=out['freewords']
vals=[mp.mpc(mp.mpf(out['vals'][1][w][0]),mp.mpf(out['vals'][1][w][1])) for w in fw]
from flint import fmpz_mat
n=110
C=mp.mpf(10)**135
rows=[]
for i,v in enumerate(vals):
    row=[0]*(n+2); row[i]=1
    row[n]=int(mp.nint(C*v.real)); row[n+1]=int(mp.nint(C*v.imag))
    rows.append(row)
B=fmpz_mat(rows).lll()
spec=[]
for i in range(B.nrows()):
    row=[int(B[i,j]) for j in range(n)]
    n2=sum(x*x for x in row)
    spec.append((n2,row))
spec.sort(key=lambda x:x[0])
def residual(row):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*vals[i]
    return abs(s)
import math
print('idx  len^2        height    log10(resid)')
for k,(n2,row) in enumerate(spec):
    r=residual(row)
    lr = -999 if r==0 else float(mp.log10(r))
    if k<35 or lr>-125:
        print(f'{k:3d} {float(n2):.3e} {max(abs(x) for x in row):.3e} {lr:9.2f}')
# count vectors with resid < 1e-131
kept=[row for n2,row in spec if residual(row)<mp.mpf('1e-131')]
print('kept (resid<1e-131):', len(kept))
import numpy as np
M=np.array(kept,dtype=object)
print('rank of kept:', np.matrix_rank(np.array(kept,dtype=float)))
