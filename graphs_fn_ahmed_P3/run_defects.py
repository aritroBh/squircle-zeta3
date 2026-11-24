"""
run_defects.py -- compute D_n = 2 - length(graph f^n) to >= 60 significant digits for the test functions,
for n = 10^3 ... 10^20 (and intermediate values for the fits). Saves data_defects.json.

Quadrature: D_n = (1/n) int_0^{vmax} g(n f(1-v/n)^{n-1} f'(1-v/n)) dv  (tanh-sinh, breakpoints around
the transition v ~ log(n a)/a), with the rigorous far-tail bound int_{x<1-vmax/n} g(y') dx <= f(1-vmax/n)^n.
Working precision dps = 70 + 2 log10(n).  Each value is computed twice (dps and dps+20) and the
difference is stored as an internal error estimate.
"""
import json, sys, time
import mpmath as mp
from multiprocessing import Pool
from fnlib import TESTS, make_callables, defect
import sympy as sp

NAMES = ['x', 'x^2', 'expm1', 'x+x^3', 'mobius', 'sin4']

def nlist():
    ns = []
    for e in range(3, 21):
        for m in [1, 2, 5]:
            ns.append(m*10**e)
    ns = [n for n in ns if n <= 10**20]
    return ns

def job(args):
    name, n = args
    expr = TESTS[name]
    f, fp = make_callables(expr)
    a_exact = sp.diff(expr, sp.Symbol('x', positive=True)).subs(sp.Symbol('x', positive=True), 1)
    dps = 70 + 2*len(str(n))
    out = []
    for extra in (0, 20):
        mp.mp.dps = dps + extra
        a = mp.mpf(sp.N(a_exact, mp.mp.dps + 10))
        t0 = time.time()
        D, bound = defect(f, fp, n, a, dps_target=dps + extra)
        out.append((D, bound, time.time() - t0))
    mp.mp.dps = dps + 20
    D0, D1 = out[0][0], out[1][0]
    return dict(name=name, n=str(n), D=mp.nstr(D1, dps + 10), est_err=mp.nstr(abs(D1 - D0), 5),
                tail_bound=mp.nstr(out[1][1], 5), secs=out[0][2] + out[1][2])

if __name__ == '__main__':
    tasks = [(name, n) for name in NAMES for n in nlist()]
    with Pool(12) as pool:
        res = pool.map(job, tasks, chunksize=1)
    json.dump(res, open('data_defects.json', 'w'), indent=1)
    for r in res:
        print(r['name'], r['n'], r['D'][:40], 'err~', r['est_err'], 'tail<=', r['tail_bound'], '%.1fs' % r['secs'])
