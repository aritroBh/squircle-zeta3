import os; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 200
r2=mp.sqrt(2); kap1=r2-1; kap2=r2/2-1; ln=mp.log
z3=mp.zeta(3); l2=ln(2); lu=ln(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
vals=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
names=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%dl2'%i for i in range(8)]+['Li2_%dlu'%i for i in range(8)]
from flint import fmpz_mat
def lll_rows(vs, dps, se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); n=len(vs)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vs[i]))
    R=fmpz_mat(rows).lll(); out=[[int(R[i,j]) for j in range(R.ncols())] for i in range(R.nrows())]
    mp.mp.dps=old; return out
# Stage 1: internal relations among basis
rows=lll_rows(vals,200,190)
rels=[r for r in rows if abs(r[-1])<100]  # residB~0 => exact relations
print("internal exact relations among basis:", len(rels))
import sympy as sp
M=sp.Matrix([r[:-1] for r in rels])
rref=M.rref()[0]
# pivot columns = dependent; free = independent
pivots=set()
for i in range(rref.rows):
    for j in range(rref.cols):
        if rref[i,j]!=0: pivots.add(j); break
indep=[i for i in range(len(vals)) if i not in pivots]
print("independent single-Li3 basis size:", len(indep), "of", len(vals))
print("independent names:", [names[i] for i in indep])
import pickle; pickle.dump({'indep':indep,'names':names}, open('indep_singles.pkl','wb'))
flog("step: reduced single-Li3 basis to %d independent"%len(indep))
