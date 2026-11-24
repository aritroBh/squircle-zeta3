import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 150
r2 = mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
J21 = mp.quad(lambda t: mp.log(1+kap1*t)**2*mp.log(1+kap2*t), [0, mp.mpf('0.5'), 1])
J12 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t)**2, [0, mp.mpf('0.5'), 1])
z3=mp.zeta(3); l2=mp.log(2); lu=mp.log(1+r2)
real_pts=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
basis=[z3]+[mp.polylog(3,p) for p in real_pts]+[mp.pi**2*l2, mp.pi**2*lu, l2**3, l2**2*lu, l2*lu**2, lu**3]
bn=['z3']+['Li3(%.4g)'%float(p) for p in real_pts]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']
from flint import fmpz_mat
def hunt(target, basis, bn, dps, scale_exp):
    old=mp.mp.dps; mp.mp.dps=dps
    B=mp.power(10,scale_exp); vals=[target]+list(basis); n=len(vals)
    rows=[]
    for i in range(n):
        row=[0]*(n+1); row[i]=1; row[-1]=int(mp.nint(B*vals[i])); rows.append(row)
    R=fmpz_mat(rows).lll()
    out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0:  # involves target
            resid=abs(v[-1]); maxc=max(abs(c) for c in v[:-1])
            out.append((resid,maxc,v))
    mp.mp.dps=old
    out.sort()
    return out
for nm,J in [('J21',J21),('J12',J12)]:
    res=hunt(J,basis,bn,150,140)
    print(nm, ": relations involving target (residB, maxc):")
    for resid,maxc,v in res[:4]:
        # show the combo
        combo=[(bn[i-1],v[i]) for i in range(1,len(v)-1) if v[i]!=0]
        print(f"   Jcoeff={v[0]} residB={resid:.2e} maxc={maxc} :: {combo}")
log("step: proper LLL hunt for J21,J12 over single-Li3 basis")
