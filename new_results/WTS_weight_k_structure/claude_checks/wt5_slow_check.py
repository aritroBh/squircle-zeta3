"""M1-session independent check of WT5 values E1, E2 (adapted from the WT4 check)'s closed forms for D1, D2 (no fastred, no TSVs, no certificates):
degree-3 kernels from the paper's SLOW sympy reducer at k = 4, 5, 6, then mpmath quadrature of
D1 = int (K5 - K4) logs and D2 = int (K6 - 2 K5 + K4) logs, compared with WT4's closed forms."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v
import mpmath as mp

K = {k: reduce_k(k, verbose=False) for k in (5, 6, 7)}
mp.mp.dps = 45
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1; lam = mp.log(2); L = mp.log(a)
tv = lambda x: (x * x + 2 * x - 1) / (4 * x)


def integral(expr_by_key):
    tot = mp.mpf(0)
    for (i, j), R in expr_by_key.items():
        f = sp.lambdify(v, R, 'mpmath')
        tot += mp.re(mp.quad(lambda x: f(x) * mp.log(tv(x)) ** i * mp.log(1 - tv(x)) ** j, [b, 1, a]))
    return tot


keys = [(4, 0), (3, 1), (2, 2)]
E1k = {kk: sp.together(K[6].get(kk, 0) - K[5].get(kk, 0)) for kk in keys}
E2k = {kk: sp.together(K[7].get(kk, 0) - 2 * K[6].get(kk, 0) + K[5].get(kk, 0)) for kk in keys}
E1n, E2n = integral(E1k), integral(E2k)
W = '<project>/research/08_theory/WT5_weight5_all_k/'
ref = lambda f: mp.mpf(open(W + f).read().split('=')[-1].strip()[:60])
print('dps', mp.mp.dps)
print('E1 slow-quad =', mp.nstr(E1n, 35), ' minus WT5 value =', mp.nstr(E1n - ref('e1_1500.txt'), 4))
print('E2 slow-quad =', mp.nstr(E2n, 35), ' minus WT5 value =', mp.nstr(E2n - ref('e2_1500.txt'), 4))
