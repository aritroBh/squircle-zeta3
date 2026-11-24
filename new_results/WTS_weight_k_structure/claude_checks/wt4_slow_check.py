"""M1-session independent check of WT4's closed forms for D1, D2 (no fastred, no TSVs, no certificates):
degree-3 kernels from the paper's SLOW sympy reducer at k = 4, 5, 6, then mpmath quadrature of
D1 = int (K5 - K4) logs and D2 = int (K6 - 2 K5 + K4) logs, compared with WT4's closed forms."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v
import mpmath as mp

K = {k: reduce_k(k, verbose=False) for k in (4, 5, 6)}
mp.mp.dps = 60
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1; lam = mp.log(2); L = mp.log(a)
tv = lambda x: (x * x + 2 * x - 1) / (4 * x)


def integral(expr_by_key):
    tot = mp.mpf(0)
    for (i, j), R in expr_by_key.items():
        f = sp.lambdify(v, R, 'mpmath')
        tot += mp.re(mp.quad(lambda x: f(x) * mp.log(tv(x)) ** i * mp.log(1 - tv(x)) ** j, [b, 1, a]))
    return tot


keys = [(3, 0), (2, 1)]
D1k = {kk: sp.together(K[5][kk] - K[4][kk]) for kk in keys}
D2k = {kk: sp.together(K[6][kk] - 2 * K[5][kk] + K[4][kk]) for kk in keys}
D1n, D2n = integral(D1k), integral(D2k)
pi2, z3, Li3 = mp.pi ** 2, mp.zeta(3), lambda z: mp.polylog(3, z)
D2cf = -lam**4/24 - pi2*lam**2/24 - mp.pi**4/160 - z3*lam/4
D1cf = (-mp.mpf(5)/48*lam**4 - lam**3*L/12 - lam**2*L**2/4 + lam*L**3/6 + L**4/3 - mp.mpf(7)/96*pi2*lam**2
        + mp.mpf(5)/24*pi2*lam*L - mp.mpf(5)/24*pi2*L**2 - mp.pi**4/192 - mp.mpf(3)/8*z3*lam - mp.mpf(7)/16*z3*L
        - lam*Li3(b) - 2*L*Li3(b) - lam*Li3(-b) + 2*L*Li3(-b) + 4*L*Li3(1/r2))
print('dps', mp.mp.dps)
print('D1 slow-quad =', mp.nstr(D1n, 40), ' closed form diff =', mp.nstr(D1n - D1cf, 5))
print('D2 slow-quad =', mp.nstr(D2n, 40), ' closed form diff =', mp.nstr(D2n - D2cf, 5))
print('D2 kernels:', {str(kk): str(sp.apart(sp.radsimp(D2k[kk]), v, extension=[sp.sqrt(2)])) for kk in keys})
