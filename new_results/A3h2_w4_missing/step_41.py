import sys, pickle, time, json
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
C=mp.mpf(10)**300
rows=[]
for i,v in enumerate(V[1]):
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
def res_exact(row,k):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V[k][i]
    return abs(s)
cands=[]
for n2,row in spec:
    if float(mp.log10(res_exact(row,1)))<-250: cands.append(row)
print('candidates:',len(cands))
# rank over Q via sympy
import sympy as sp
M=sp.Matrix(cands)
print('rank of 29 candidates:', M.rank())
# residuals at all 4 embeddings (classify form-level vs exotic)
print('idx support height  log10res k=1   k=3   k=5   k=7')
for i,row in enumerate(cands):
    rs=[res_exact(row,k) for k in (1,3,5,7)]
    supp=sum(1 for x in row if x)
    print(f'{i:3d} {supp:4d} {max(abs(x) for x in row):7.1e} '+' '.join(f'{-999 if r==0 else float(mp.log10(r)):8.1f}' for r in rs))
pickle.dump({'freewords':fw,'cands':cands},open(WORK+'/cands29.pickle','wb'))
print('saved cands29.pickle')
