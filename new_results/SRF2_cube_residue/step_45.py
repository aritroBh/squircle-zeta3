import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import ts_nodes
from srf2_lib import Phi, tomp
ctx.dps = 80; mp.mp.dps = 80

def E_int(eps, h=arb(1)/32, T=4.6):
    """E(eps) = -(1+2e)/(2e+1) + int_0^1 a^{2e-1}[(1+a)^{-1-2e} S0 - 1 + (1+2e) a] da."""
    q = 1+2*eps
    nodes = ts_nodes(h, T)
    s = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = Phi(a, q)*S0 + (1-q*a)*A2/(S0+1)
        s += (arb(2*eps)*la).exp()/a*f*wa
    return -arb(1+2*eps)/(2*eps+1) + s

for e in ['-0.9', '-0.95', '-0.99', '-1', '-0.999']:
    v = E_int(arb(e))
    print(f"E({e}) = {mp.nstr(tomp(v), 25)}")
mp.mp.dps = 60
E_m1 = -1 + mp.quad(lambda a: ((1+a)*mp.sqrt(1+a**4)-1-a)/a**3, [0, mp.mpf('0.5'), 1])
print("E(-1) via mpmath =", mp.nstr(E_m1, 25))
# cross-check E at eps=0.1 vs vlib E_series
from vlib import E_series, tompf
Ef, _ = E_series(2, 0.1, 2, h=arb(1)/32, T=4.6)
print("raw E(0.1) =", mp.nstr(tompf(Ef[0],90), 20), "  E_int(0.1) =", mp.nstr(tomp(E_int(arb('0.1'))), 20))
