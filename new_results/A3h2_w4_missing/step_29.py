import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=200
BASE='<project>'
vh=pickle.load(open(BASE+'/research/08_theory/A3f4_glanois_basis/vals_hi.pickle','rb'))
five=pickle.load(open(BASE+'/research/08_theory/A3f4_glanois_basis/five_relations.pickle','rb'))
fw=five['freewords']; R5=five['relations']
vals=[mp.mpc(mp.mpf(vh['vals'][27+i][0]), mp.mpf(vh['vals'][27+i][1])) for i in range(32)]
from flint import fmpz_mat
C=mp.mpf(10)**170
n=32
rows=[]
for i,v in enumerate(vals):
    row=[0]*(n+2); row[i]=1
    row[n]=int(mp.nint(C*v.real)); row[n+1]=int(mp.nint(C*v.imag))
    rows.append(row)
L=fmpz_mat(rows)
t=time.time(); B=L.lll(); print('flint lll',time.time()-t)
norms=[]
for i in range(B.nrows()):
    row=[int(B[i,j]) for j in range(32)]
    tail=[int(B[i,j]) for j in (32,33)]
    n2=sum(x*x for x in row)
    norms.append((n2,row,tail))
norms.sort(key=lambda x:x[0])
def residual(row):
    s=mp.mpc(0)
    for i,c in enumerate(row): s+=c*vals[i]
    return abs(s)
for n2,row,tail in norms[:9]:
    print('len^2=%.3e height=%d resid=%s tail=%s'%(float(n2),max(abs(x) for x in row),mp.nstr(residual(row),3),[mp.nstr(mp.mpf(t)/C,3) for t in tail]))
print('---R1..R5---')
for i,r in enumerate(R5):
    print(f'R{i+1}:', mp.nstr(residual(r),3), 'height', max(abs(x) for x in r))
