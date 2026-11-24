import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 150
r2 = mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
# cleanest weight-3 cross integral
K = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t)/t, [mp.mpf(0), mp.mpf('0.5'), 1])
print("K =", mp.nstr(K,30))
# comprehensive single-polylog weight-3 basis
z3=mp.zeta(3); l2=mp.log(2); lu=mp.log(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]
Li2=[mp.polylog(2,p) for p in reps]
basis=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]
# add Li2*log products
basis += [Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
bn=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%d l2'%i for i in range(8)]+['Li2_%d lu'%i for i in range(8)]
from flint import fmpz_mat
def hunt(target, basis, dps, scale_exp):
    old=mp.mp.dps; mp.mp.dps=dps
    B=mp.power(10,scale_exp); vals=[target]+list(basis); n=len(vals)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n):
        rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vals[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0:
            out.append((abs(v[-1]), max(abs(c) for c in v[:-1]), v))
    mp.mp.dps=old; out.sort(); return out
res=hunt(K,basis,150,135)
print("relations involving K (residB, maxc):")
for resid,maxc,v in res[:5]:
    combo=[(bn[i-1],v[i]) for i in range(1,len(v)-1) if v[i]!=0]
    print(f"  Kc={v[0]} residB={resid:.1e} maxc={maxc} :: {combo[:6]}")
log("step: tested K against fuller single-polylog basis")
