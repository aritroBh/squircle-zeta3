"""Claude's independent check for WT3: weight-3 part of M_k from the SLOW sympy reducer (paper's ibp_reduce.reduce_k),
evaluated by direct mpmath quadrature. Shares no code with WT2/WT3's fast reducer.
Usage: python wt3_slow_check.py k1 k2 ...   (writes w3_slow_<k>.json next to this script)
"""
import sys, json, os, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 50
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1; lam = mp.log(2); L = mp.log(a)


def tv(x):
    return (x * x + 2 * x - 1) / (4 * x)


def w3_known(k):
    z3 = mp.zeta(3); pi2 = mp.pi ** 2
    W3p = (7*r2/32 - mp.mpf(3)/16)*z3 + (1+r2)*mp.polylog(3, b) + (1-r2)*mp.polylog(3, -b) \
        - 2*r2*mp.polylog(3, 1/r2) + (2*L - r2*lam)*mp.polylog(2, b) \
        + (r2/12 - mp.mpf(1)/8)*pi2*L - (1+r2)/24*pi2*lam + (1+r2)/24*lam**3 \
        + r2/4*lam**2*L - (1+r2)/2*lam*L**2 + (5-r2)/6*L**3
    if k == 3:
        return W3p
    if k == 4:
        return W3p - z3/8 - pi2*lam/24 - lam**3/12
    return None


for k in [int(x) for x in sys.argv[1:]]:
    t0 = time.time()
    out = reduce_k(k, verbose=False)
    keys = sorted(kk for kk in out if sum(kk) == 2)
    total = mp.mpf(0)
    parts = {}
    for (i, j) in keys:
        f = sp.lambdify(v, out[(i, j)], 'mpmath')
        def integrand(x, f=f, i=i, j=j):
            t = tv(x)
            return f(x) * mp.log(t) ** i * mp.log(1 - t) ** j
        val = mp.quad(integrand, [b, 1, a])
        val = mp.re(val)
        parts[f'{i},{j}'] = mp.nstr(val, 45)
        total += val
    known = w3_known(k)
    res = {'k': k, 'dps': mp.mp.dps, 'keys': [list(kk) for kk in keys], 'parts': parts,
           'W3_slow': mp.nstr(total, 45),
           'kernels': {f'{i},{j}': str(out[(i, j)]) for (i, j) in keys},
           'seconds': round(time.time() - t0, 1)}
    if known is not None:
        res['diff_vs_paper'] = mp.nstr(total - known, 5)
    print(json.dumps({x: res[x] for x in res if x != 'kernels'}), flush=True)
    with open(os.path.join(HERE, f'w3_slow_{k}.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
