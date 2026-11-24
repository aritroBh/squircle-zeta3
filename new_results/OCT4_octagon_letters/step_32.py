import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 120
r2 = mp.sqrt(2)
kap1 = r2-1; kap2 = r2/2-1
print("kap1=",mp.nstr(kap1,10)," kap2=",mp.nstr(kap2,10))
# weight-3 double-log integrals
J21 = mp.quad(lambda t: mp.log(1+kap1*t)**2*mp.log(1+kap2*t), [0, mp.mpf('0.5'), 1])
J12 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t)**2, [0, mp.mpf('0.5'), 1])
print("J21 =", mp.nstr(J21,30))
print("J12 =", mp.nstr(J12,30))
# single-Li3 weight-3 basis
z3=mp.zeta(3); l2=mp.log(2); lu=mp.log(1+r2)
real_pts=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
basis=[z3]+[mp.polylog(3,p) for p in real_pts]+[mp.pi**2*l2, mp.pi**2*lu, l2**3, l2**2*lu, l2*lu**2, lu**3]
bn=['z3']+['Li3(%s)'%p for p in real_pts]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']
# hunt J21, J12 over this basis (rational coeffs) via pslq
for nm,J in [('J21',J21),('J12',J12)]:
    rel = mp.pslq([J]+basis, tol=mp.mpf('1e-100'), maxcoeff=10**8, maxsteps=10000)
    if rel is None:
        print(nm, ": no relation over single-Li3 basis (maxcoeff 1e8)")
    else:
        print(nm, "relation:", rel[0], [ (bn[i-1],rel[i]) for i in range(1,len(rel)) if rel[i]])
log("step: J21,J12 double-log integrals do NOT reduce to single Li3's (pslq negative)")
