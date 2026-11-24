import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import W_int, J_int, tomp
ctx.dps = 60; mp.mp.dps = 60

# J(-1) exactly: J1(-1) = int a^{-2}(sqrt(1+a^4)-1) da (smooth), J(-1) = -1 + J1(-1)
Jm1 = J_int(arb(-1))
print("J(-1) =", mp.nstr(tomp(Jm1), 25))
I_closed = 1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("-1 + I =", mp.nstr(-1 + I_closed, 25))

# h-convergence of J near eps=-1
for e in ['-0.99', '-0.999']:
    for hh in ['32','64']:
        eps = arb(e)
        J = J_int(eps, h=arb(1)/int(hh), T=5.2)
        print(f"J({e}) h=1/{hh}: {mp.nstr(tomp(J), 20)}")
