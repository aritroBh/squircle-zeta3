import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=160
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
t=time.time(); B=fmpz_mat(rows).lll(); print('flint 8-constraint lll:',time.time()-t,'s')
spec=[]
for i in range(B.nrows()):
    row=[int(B[i,j]) for j in range(n)]
    n2=sum(x*x for x in row)
    spec.append((n2,row))
spec.sort(key=lambda x:x[0])
def residual_all(row):
    rs=[]
    for k in (1,3,5,7):
        s=mp.mpc(0)
        for i,c in enumerate(row):
            if c: s+=c*V[k][i]
        rs.append(abs(s))
    return rs
print('idx len^2     height  log10resid(k=1,3,5,7)')
for k,(n2,row) in enumerate(spec[:40]):
    rs=residual_all(row)
    print(f'{k:3d} {float(n2):9.2e} {max(abs(x) for x in row):8.2e} '+' '.join(f'{-999 if r==0 else float(mp.log10(r)):8.1f}' for r in rs))
