import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=340
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
fw=out['freewords']
V={k:[mp.mpc(mp.mpf(out['vals'][k][w][0]),mp.mpf(out['vals'][k][w][1])) for w in fw] for k in (1,3,5,7)}
from flint import fmpz_mat
n=110
C=mp.mpf(10)**105
rows=[]
for i in range(n):
    row=[0]*(n+8); row[i]=1
    j=n
    for k in (1,3,5,7):
        row[j]=int(mp.nint(C*V[k][i].real)); row[j+1]=int(mp.nint(C*V[k][i].imag)); j+=2
    rows.append(row)
B=fmpz_mat(rows).lll()
spec=[]
for i in range(B.nrows()):
    row=[int(B[i,j]) for j in range(n)]
    n2=sum(x*x for x in row)
    spec.append((n2,row))
spec.sort(key=lambda x:x[0])
def res_exact(row,k=1):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V[k][i]
    return abs(s)
print('idx height   log10 exact-resid k=1 (true<-250)')
cnt=0
for kk,(n2,row) in enumerate(spec):
    r=res_exact(row,1)
    lr = -999 if r==0 else float(mp.log10(r))
    if lr < -250: cnt+=1
    if kk<45:
        print(f'{kk:3d} {max(abs(x) for x in row):9.2e} {lr:10.1f}')
print('count with exact-resid<1e-250:', cnt)
