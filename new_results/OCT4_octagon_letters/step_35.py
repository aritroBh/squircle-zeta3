import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 150
r2 = mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
def Jj(pol, a, b):
    return mp.quad(lambda t: pol(t)*mp.log(1+kap1*t)**a*mp.log(1+kap2*t)**b, [mp.mpf(0),mp.mpf('0.5'),1])
one=lambda t: 1; tt=lambda t:t; t2=lambda t:t*t; omt=lambda t:1-t
J21=Jj(one,2,1); J12=Jj(one,1,2)
J21_t=Jj(tt,2,1); J12_t=Jj(tt,1,2); J21_1mt=Jj(omt,2,1); J12_1mt=Jj(omt,1,2)
z3=mp.zeta(3); l2=mp.log(2); lu=mp.log(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
basis=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
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
        if v[0]!=0: out.append((abs(v[-1]), max(abs(c) for c in v[:-1]), v))
    mp.mp.dps=old; out.sort(); return out
for nm,J in [('J21',J21),('J12',J12),('J21_t',J21_t),('J12_t',J12_t),('J21_1mt',J21_1mt),('J12_1mt',J12_1mt)]:
    res=hunt(J,basis,150,135)
    good=[r for r in res if r[1]<10000]
    if good:
        r=good[0]; print(f"{nm}: SMALL relation! maxc={r[1]} residB={r[0]:.1e}")
    else:
        print(f"{nm}: only spurious (floor maxc={res[0][1]:.1e})" if res else f"{nm}: none")
log("step: re-hunted cross integrals with fuller basis")
