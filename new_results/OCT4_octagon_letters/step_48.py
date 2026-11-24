import os, pickle; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 280
d=pickle.load(open('basis280.pkl','rb'))
allv=[mp.mpf(x) for x in d['allv']]; alln=d['alln']
from flint import fmpz_mat
def lll_rows(vs,dps,se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); n=len(vs)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vs[i]))
    R=fmpz_mat(rows).lll(); out=[[int(R[i,j]) for j in range(R.ncols())] for i in range(R.nrows())]
    mp.mp.dps=old; return out
rows=lll_rows(allv,280,265)
rels=[r for r in rows if abs(r[-1])<1000]
print("internal exact relations:", len(rels))
import sympy as sp
rref=sp.Matrix([r[:-1] for r in rels]).rref()[0]
pivots=set()
for i in range(rref.rows):
    for j in range(rref.cols):
        if rref[i,j]!=0: pivots.add(j); break
indep=[i for i in range(len(allv)) if i not in pivots]
print("independent basis size:", len(indep), "of", len(allv))
print("indep names:", [alln[i] for i in indep])
pickle.dump(indep, open('indep280.pkl','wb'))
flog("step: reduced 46 consts to %d independent at 280dps"%len(indep))
