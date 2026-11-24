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
    tail=[int(B[i,j]) for j in (n,n+1)]
    n2=sum(x*x for x in row)
    t2=sum(x*x for x in tail)
    spec.append((n2,t2,row,tail))
spec.sort(key=lambda x:x[0]+x[1])
print('total rows:', B.nrows())
print('idx   log10|a|      log10|tail|   log10|tail/C|')
for k,(n2,t2,row,tail) in enumerate(spec):
    la = 0.5*float(mp.log10(max(n2,1)))
    lt = 0.5*float(mp.log10(max(t2,1)))
    print(f'{k:3d} {la:10.2f} {lt:10.2f} {lt-135:10.2f}')
