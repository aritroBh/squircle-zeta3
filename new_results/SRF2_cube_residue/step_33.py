import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import J2_int, W2_int, tomp
ctx.dps = 60; mp.mp.dps = 60

# h-convergence of J2 at -0.999
for hh, TT in [('32',4.6), ('64',5.2), ('128',5.6)]:
    v = J2_int(arb('-0.999'), h=arb(1)/int(hh), T=TT)
    print(f"J2(-0.999) h=1/{hh} T={TT}: {mp.nstr(tomp(v), 22)}")

# independent mpmath quad of J2's integral part at eps=-0.999
mp.mp.dps = 40
eps = mp.mpf('-0.999')
def integrand(a):
    if a == 0: return mp.mpf('0.001002')*0  # a^{-1.998}*O(a^2) -> 0... handle below
    S0 = mp.sqrt(1+a**(2-2*eps))
    f1 = (1+a)**(-2-2*eps)*S0 - 1
    return a**(2*eps)*(f1 + (2+2*eps)*a)
# split [0,1] to help quad
J2int = mp.quad(integrand, [0, mp.mpf('1e-3'), mp.mpf('0.01'), mp.mpf('0.1'), mp.mpf('0.5'), 1])
print("mpmath J2_int part:", mp.nstr(J2int, 22))
J2_full = 1/(2*eps+1) - (2+2*eps)/(2*eps+2) + J2int
print("mpmath J2(-0.999) =", mp.nstr(J2_full, 22))
print("J(-1) value       = -0.84721308479397908661")
