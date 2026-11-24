import os, pickle; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 200
r2=mp.sqrt(2); kap1=r2-1; kap2=r2/2-1; ln=mp.log
x1=1-r2; x2=1-r2/2
def Li21(x,y):
    x=mp.mpf(x); y=mp.mpf(y); N=int(mp.mp.dps/max(mp.log10(1/abs(x)),mp.mpf('0.05')))+60
    s=mp.mpf(0); Hy=mp.mpf(0); py=mp.mpf(1)
    for n in range(1,N+1):
        s += x**n/n**2*Hy; py*=y; Hy+=py/n
    return s
cand={'Li21(x1,x2)':(x1,x2),'Li21(x2,x1)':(x2,x1),'Li21(x1,x2/x1)':(x1,x2/x1),'Li21(x2/x1,x1)':(x2/x1,x1),
      'Li21(x2,x2/x1)':(x2,x2/x1),'Li21(x2/x1,x2)':(x2/x1,x2)}
z3=mp.zeta(3); l2=ln(2); lu=ln(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
vals=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
names=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%dl2'%i for i in range(8)]+['Li2_%dlu'%i for i in range(8)]
indep=pickle.load(open('indep_singles.pkl','rb'))['indep']
bvals=[vals[i] for i in indep]; bn=[names[i] for i in indep]
T_t  = mp.quad(lambda t: ln(1+kap1*t)*ln(1+kap2*t)*ln(t), [mp.mpf(0),mp.mpf('0.25'),mp.mpf('0.75'),1])
from flint import fmpz_mat
def hunt1(target,bvals,extra,dps,se):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se); vs=[target]+list(bvals)+[extra]; n=len(vs)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vs[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0: out.append((abs(v[-1]),max(abs(c) for c in v[:-1]),v))
    mp.mp.dps=old; out.sort(); return out
print("baseline (no Li21):")
res=hunt1(T_t,bvals,mp.mpf('0.12345'),200,180); print("  floor maxc=%.1e"%res[0][1])
for k,(a,b) in cand.items():
    c=Li21(a,b); res=hunt1(T_t,bvals,c,200,180); r0=res[0]
    tag="CLEAN!" if r0[0]<1000 else ""
    print(f"  +{k}: residB={r0[0]:.1e} maxc={r0[1]:.1e} {tag}")
flog("step: one-at-a-time Li21 hunt for T_t")
