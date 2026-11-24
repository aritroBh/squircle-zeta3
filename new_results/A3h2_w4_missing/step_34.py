import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=180
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
fw=out['freewords']
vals=[mp.mpc(mp.mpf(out['vals'][1][w][0]),mp.mpf(out['vals'][1][w][1])) for w in fw]
print('max |val|:', mp.nstr(max(abs(v) for v in vals),4), 'min nonzero:', mp.nstr(min(abs(v) for v in vals),4))
from flint import fmpz_mat
n=110
C=mp.mpf(10)**135
rows=[]
for i,v in enumerate(vals):
    row=[0]*(n+2); row[i]=1
    row[n]=int(mp.nint(C*v.real)); row[n+1]=int(mp.nint(C*v.imag))
    rows.append(row)
t=time.time(); L=fmpz_mat(rows); print('build',time.time()-t)
t=time.time(); B=L.lll(); print('flint lll time',time.time()-t)
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
print('--- 40 shortest ---')
for n2,row in spec[:40]:
    print('len^2=%.3e height=%.3e resid=%s'%(float(n2),max(abs(x) for x in row),mp.nstr(residual(row),3)))
