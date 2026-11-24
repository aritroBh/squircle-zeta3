import os; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 120
r2=mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
ln=mp.log
# three-letter cross integrals (pure triple-log, weight 3)
T_t  = mp.quad(lambda t: ln(1+kap1*t)*ln(1+kap2*t)*ln(t), [mp.mpf(0),mp.mpf('0.25'),mp.mpf('0.75'),1])
T_1mt= mp.quad(lambda t: ln(1+kap1*t)*ln(1+kap2*t)*ln(1-t), [mp.mpf(0),mp.mpf('0.25'),mp.mpf('0.75'),1])
print("T_t  =", mp.nstr(T_t,25))
print("T_1mt=", mp.nstr(T_1mt,25))
# single-Li3 independent basis (reduce first)
z3=mp.zeta(3); l2=ln(2); lu=ln(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
basis=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
from flint import fmpz_mat
def huntJ(target,basis,dps,se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); vals=[target]+list(basis); n=len(vals)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vals[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0: out.append((abs(v[-1]),max(abs(c) for c in v[:-1]),v[0]))
    mp.mp.dps=old; out.sort(); return out
for nm,T in [('T_t',T_t),('T_1mt',T_1mt)]:
    res=huntJ(T,basis,120,112)
    # smallest-residual relation with target-coeff !=0
    print(nm,": best target-involving relations (residB, maxc, Jc):", [(r[0],r[1],r[2]) for r in res[:3]])
flog("step: tested three-letter cross integrals T_t, T_1mt for Li_2,1 need")
