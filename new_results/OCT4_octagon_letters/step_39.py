import os; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 250
r2=mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
J21 = mp.quad(lambda t: mp.log(1+kap1*t)**2*mp.log(1+kap2*t), [mp.mpf(0),mp.mpf('0.5'),1])
J12 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t)**2, [mp.mpf(0),mp.mpf('0.5'),1])
g1=r2-1; g2=2-r2  # the two predicted Li3 points
l2=mp.log(2); lu=mp.log(1+r2); z3=mp.zeta(3); pi2=mp.pi**2
Li3g1=mp.polylog(3,g1); Li3g2=mp.polylog(3,g2); Li2g1=mp.polylog(2,g1); Li2g2=mp.polylog(2,g2)
basis=[Li3g1,Li3g2,Li2g1*l2,Li2g1*lu,Li2g2*l2,Li2g2*lu,z3,pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3,pi2,l2**2,l2*lu,lu**2,Li2g1,Li2g2,l2,lu,mp.mpf(1)]
bn=['Li3g1','Li3g2','Li2g1l2','Li2g1lu','Li2g2l2','Li2g2lu','z3','pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3','pi2','l2^2','l2lu','lu^2','Li2g1','Li2g2','l2','lu','1']
for nm,J in [('J21',J21),('J12',J12)]:
    rel=mp.pslq([J]+basis, tol=mp.mpf('1e-235'), maxcoeff=10**12, maxsteps=100000)
    if rel is None: print(nm,": NO relation (maxcoeff 1e12)")
    else:
        # verify
        s=rel[0]*J+sum(rel[i+1]*basis[i] for i in range(len(basis)))
        print(nm,": rel found, Jc=%d, residual=%s, maxc=%d"%(rel[0], mp.nstr(abs(s),3), max(abs(c) for c in rel)))
        print("   ", [(bn[i-1],rel[i]) for i in range(1,len(rel)) if rel[i]!=0])
log("step: focused pslq for J21/J12")
