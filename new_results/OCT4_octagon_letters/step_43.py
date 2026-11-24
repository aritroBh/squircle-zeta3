import os, pickle; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 200
r2=mp.sqrt(2); kap1=r2-1; kap2=r2/2-1; ln=mp.log
z3=mp.zeta(3); l2=ln(2); lu=ln(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
vals=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
names=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%dl2'%i for i in range(8)]+['Li2_%dlu'%i for i in range(8)]
indep=pickle.load(open('indep_singles.pkl','rb'))['indep']
bvals=[vals[i] for i in indep]; bn=[names[i] for i in indep]
T_t  = mp.quad(lambda t: ln(1+kap1*t)*ln(1+kap2*t)*ln(t), [mp.mpf(0),mp.mpf('0.25'),mp.mpf('0.75'),1])
T_1mt= mp.quad(lambda t: ln(1+kap1*t)*ln(1+kap2*t)*ln(1-t), [mp.mpf(0),mp.mpf('0.25'),mp.mpf('0.75'),1])
from flint import fmpz_mat
def huntJ(target,bvals,dps,se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); vs=[target]+list(bvals); n=len(vs)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vs[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0: out.append((abs(v[-1]),max(abs(c) for c in v[:-1]),v))
    mp.mp.dps=old; out.sort(); return out
for nm,T in [('T_t',T_t),('T_1mt',T_1mt)]:
    res=huntJ(T,bvals,200,185)
    r0=res[0]
    print(nm,": best target relation residB=%d maxc=%d Jc=%d"%(r0[0],r0[1],r0[2][0]))
    if r0[0]<1000:
        combo=[(bn[i-1],r0[2][i]) for i in range(1,len(r0[2])-1) if r0[2][i]!=0]
        print("   RELATION:", combo)
    else:
        print("   (spurious; no exact relation)")
flog("step: hunted T_t,T_1mt vs independent singles")
