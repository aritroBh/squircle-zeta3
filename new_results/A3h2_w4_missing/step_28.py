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
import numpy as np
def build_lattice(vals2, C):
    n=len(vals2)
    rows=[]
    for i,v in enumerate(vals2):
        row=[0]*(n+2)
        row[i]=1
        row[n]=int(mp.nint(C*v.real))
        row[n+1]=int(mp.nint(C*v.imag))
        rows.append(row)
    return fmpz_mat(rows)

C=mp.mpf(10)**170
t=time.time()
L=build_lattice(vals,C)
print('build',time.time()-t)
t=time.time()
B=L.lll()
print('flint lll time',time.time()-t)
# inspect shortest vectors
norms=[]
for i in range(B.nrows()):
    row=[B[i,j] for j in range(32)]
    tail=[B[i,j] for j in (32,33)]
    n2=sum(x*x for x in row)
    norms.append((n2,row,tail))
norms.sort(key=lambda x:x[0])
for n2,row,tail in norms[:8]:
    print('len^2=%.3e'%float(n2),'height=%.3e'%max(abs(x) for x in row),'tail=',[mp.nstr(mp.mpf(t)/C,3) for t in tail])
# check R1-R5 membership
def residual(row):
    s=mp.mpc(0)
    for i,c in enumerate(row): s+=c*vals[i]
    return abs(s)
for i,r in enumerate(R5):
    print(f'R{i+1} residual:', mp.nstr(residual(r),3), 'height', max(abs(x) for x in r))
