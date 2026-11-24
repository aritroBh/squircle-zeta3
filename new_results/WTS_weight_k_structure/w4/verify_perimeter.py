"""verify_perimeter.py <dps> <n>
Sanity check of the DEFINITION chain squircle -> M_k:  the perimeter of |x|^p+|y|^p=1 is P(p) = 4 eps int_0^1 sqrt(t^{2eps-2}+(1-t)^{2eps-2}) dt, eps = 1/p
(arclength of x=t^{eps}, y=(1-t)^{eps} in the variable t = x^p).  Compare P(p) computed by direct quadrature with 8(1 + sum_{k=0}^{K} M_k p^{-(k+1)}),
M_k = my direct integrals (verify_lib.mk_direct_integrand).   PRECISION: mp.mp.dps set after all imports and printed."""
import sys
sys.dont_write_bytecode = True
from verify_lib import *
dps, n = int(sys.argv[1]), int(sys.argv[2])
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps, flush=True)
X = default_X(dps)
K = 16
M = [integrate([mk_direct_integrand(k)], X, n, halves=(0,))[0] for k in range(K + 1)]
for k in range(K + 1):
    print(f"  M_{k} = {mp.nstr(M[k], 25)}")
def P(p):
    """arclength of the first-quadrant arc, times 4:  P = 4 eps int_0^1 sqrt(a^2+b^2) dt, a = t^{eps-1}, b = (1-t)^{eps-1}.
    int_0^{1/2} sqrt(a^2+b^2) = int_0^{1/2} a dt (= (1/2)^eps/eps exactly) + int_0^{1/2} b^2/(sqrt(a^2+b^2)+a) dt  (no cancellation)."""
    eps = mp.mpf(1) / p
    def g(t):
        a = t ** (eps - 1); b = (1 - t) ** (eps - 1)
        return b * b / (mp.sqrt(a * a + b * b) + a)
    half = mp.mpf(2) ** (-eps) / eps + mp.quad(g, [0, mp.mpf(1) / 2])
    return 4 * eps * 2 * half
for p in (50, 100):
    num = P(p)
    ser = 8 * (1 + sum(M[k] * mp.mpf(p) ** (-(k + 1)) for k in range(K + 1)))
    ser0 = 8 * (1 + sum(M[k] * mp.mpf(p) ** (-(k + 1)) for k in range(K)))
    print(f"p={p}: P(p) quadrature = {mp.nstr(num, 40)};  series(M_0..M_{K}) - P = {mp.nstr(ser - num, 5)};  series(M_0..M_{K-1}) - P = {mp.nstr(ser0 - num, 5)}  (mp.dps {mp.mp.dps})")
