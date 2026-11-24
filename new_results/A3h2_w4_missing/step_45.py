import sys, pickle, random
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=680
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
d=pickle.load(open(WORK+'/cands29.pickle','rb'))
fw=d['freewords']; true29=d['cands']
V1=[mp.mpc(mp.mpf(a),mp.mpf(b)) for a,b in pickle.load(open(WORK+'/vals110_k1_700.pickle','rb'))['V1_700']]
def res(row):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V1[i]
    return abs(s)
random.seed(7)
print('--- proper lattice combos (should vanish ~1e-676) ---')
for trial in range(3):
    co=[random.randint(-50,50) for _ in range(29)]
    row=[sum(co[j]*true29[j][i] for j in range(29)) for i in range(110)]
    r=res(row)
    print(f'combo {trial}: log10res={float(mp.log10(r)) if r>0 else -999:8.1f}')
# rank of the 29 mod p
from flint import nmod_mat
p=1000000009
M=nmod_mat(29,110,p)
for i,row in enumerate(true29):
    for j,x in enumerate(row):
        if x: M[i,j]=x%p
print('rank of 29 relations mod p:', M.rank())
# LLL-reduce the relation lattice to a short basis
from flint import fmpz_mat
B=fmpz_mat([list(r) for r in true29]).lll()
red=[[int(B[i,j]) for j in range(110)] for i in range(B.nrows())]
print('reduced basis size:',len(red))
print('--- reduced basis: heights, len^2, support, resid ---')
for i,row in enumerate(red):
    r=res(row)
    print(f'{i:3d} height={max(abs(x) for x in row):8.2e} len2={sum(x*x for x in row):10.2e} supp={sum(1 for x in row if x):4d} log10res={float(mp.log10(r)) if r>0 else -999:7.1f}')
pickle.dump({'freewords':fw,'cands':red},open(WORK+'/cands29_reduced.pickle','wb'))
