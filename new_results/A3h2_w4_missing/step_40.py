import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=340
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
fw=out['freewords']
V1=[mp.mpc(mp.mpf(out['vals'][1][w][0]),mp.mpf(out['vals'][1][w][1])) for w in fw]
from flint import fmpz_mat
n=110
C=mp.mpf(10)**300
rows=[]
for i,v in enumerate(V1):
    row=[0]*(n+2); row[i]=1
    row[n]=int(mp.nint(C*v.real)); row[n+1]=int(mp.nint(C*v.imag))
    rows.append(row)
t=time.time(); B=fmpz_mat(rows).lll(); print('flint lll C=1e300:',f'{time.time()-t:.1f}s')
spec=[]
for i in range(B.nrows()):
    row=[int(B[i,j]) for j in range(n)]
    n2=sum(x*x for x in row)
    spec.append((n2,row))
spec.sort(key=lambda x:x[0])
def res_exact(row):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V1[i]
    return abs(s)
cands=[]
print('idx height   log10 exact-resid')
for kk,(n2,row) in enumerate(spec):
    r=res_exact(row)
    lr=-999 if r==0 else float(mp.log10(r))
    if lr<-250: cands.append(row)
    if kk<40 or lr<-250: print(f'{kk:3d} {max(abs(x) for x in row):9.2e} {lr:10.1f}')
print('TRUE relations found:',len(cands))
