import os, sys, json, pickle; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 120
r2 = mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
cv = {k: mp.mpf(v) for k,v in pickle.load(open('li21_cand.pkl','rb')).items()}
def Jj(pol,a,b): return mp.quad(lambda t: pol(t)*mp.log(1+kap1*t)**a*mp.log(1+kap2*t)**b,[mp.mpf(0),mp.mpf('0.5'),1])
one=lambda t:1; J21=Jj(one,2,1); J12=Jj(one,1,2)
z3=mp.zeta(3); l2=mp.log(2); lu=mp.log(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
singles=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
sn=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%d.l2'%i for i in range(8)]+['Li2_%d.lu'%i for i in range(8)]
basis=singles+[cv[k] for k in cv]; bn=sn+list(cv.keys())
from flint import fmpz_mat
def hunt(target,basis,dps,se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); vals=[target]+list(basis); n=len(vals)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vals[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0: out.append((abs(v[-1]),max(abs(c) for c in v[:-1]),v))
    mp.mp.dps=old; out.sort(); return out
for nm,J in [('J21',J21),('J12',J12)]:
    res=hunt(J,basis,120,110); good=[r for r in res if r[1]<10**6]
    print(nm,":", ("SMALL relations:" if good else "only spurious floor maxc=%.1e"%(res[0][1] if res else 0)))
    for r in good[:3]:
        combo=[(bn[i-1],r[2][i]) for i in range(1,len(r[2])-1) if r[2][i]!=0]
        print("   maxc=%d residB=%.1e ::"%(r[1],r[0]), combo)
log("step: hunted J21/J12 with Li21 candidates added")
