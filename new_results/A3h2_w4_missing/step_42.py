import sys, pickle, time, json
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
t=time.time(); B=fmpz_mat(rows).lll(); print('lll',f'{time.time()-t:.1f}s',flush=True)
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
t=time.time()
cands=[]
for n2,row in spec:
    r=res_exact(row)
    if r>0 and float(mp.log10(r))<-250: cands.append(row)
    elif r==0: cands.append(row)
print('resid loop',f'{time.time()-t:.1f}s','cands',len(cands),flush=True)
pickle.dump({'freewords':fw,'cands':cands},open(WORK+'/cands29.pickle','wb'))
print('saved')
