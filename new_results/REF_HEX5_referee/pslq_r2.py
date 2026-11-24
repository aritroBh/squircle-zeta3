"""Independent check of the r = 2 example (paper hexagon |x|^p+|y|^p+3|x-y|^p <= 4^p, weighted convention):
integer relation for X = E_e/sqrt5 of the r = 2 edge (own direct jets) on the 20 alternating-MZV constants.
Compare its zeta(5), zeta(2)zeta(3) coefficients with the paper's -20187/2000 and -858/125 (and the
local law -2079/250, -477/250), and check that Li4(1/2), Li5(1/2) occur."""
import sys, time, mpmath as mp, sympy as sp
import direct5 as D, formula5 as FM
from run_hex5 import CASES
names = ['1', 'l2', 'pi2', 'l2^2', 'z3', 'pi2 l2', 'l2^3', 'pi4', 'z3 l2', 'pi2 l2^2', 'l2^4', 'Li4h',
         'z5', 'z2z3', 'pi4 l2', 'z3 l2^2', 'pi2 l2^3', 'l2^5', 'Li5h', 'Li4h l2']
res = {}
for dps in [int(a) for a in sys.argv[1:]]:
    mp.mp.dps = dps
    t0 = time.time()
    V, w = CASES['paper_2_half_1']
    H = D.setup([(FM.N(a), FM.N(b)) for a, b in V], [FM.N(x) for x in w])
    E, err, lev = D.edge_derivs(H, 0)
    X = E[5] / mp.sqrt(5)
    l2 = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3); z5 = mp.zeta(5)
    L4 = mp.polylog(4, mp.mpf(1) / 2); L5 = mp.polylog(5, mp.mpf(1) / 2)
    basis = [1, l2, pi2, l2 ** 2, z3, pi2 * l2, l2 ** 3, pi2 ** 2, z3 * l2, pi2 * l2 ** 2, l2 ** 4, L4,
             z5, pi2 / 6 * z3, pi2 ** 2 * l2, z3 * l2 ** 2, pi2 * l2 ** 3, l2 ** 5, L5, L4 * l2]
    print('dps %d: E = %s (err est %s, %.0f s)' % (dps, mp.nstr(E[5], 40), mp.nstr(err, 3), time.time() - t0))
    rel = mp.pslq([X] + basis, maxcoeff=10 ** 13, maxsteps=10 ** 7, tol=mp.mpf(10) ** (-(dps - 20)))
    if rel is None:
        print('  no relation'); continue
    resid = mp.fsum(c * v for c, v in zip(rel, [X] + basis))
    c = {nm: sp.Rational(-rel[i + 1], rel[0]) for i, nm in enumerate(names)}
    print('  relation maxcoeff %d, residual %s (%.0f s)' % (max(abs(x) for x in rel), mp.nstr(resid, 3), time.time() - t0))
    print('  coefficients of X:', c)
    print('  [z5] = %s (paper -20187/2000) ; [z2z3] = %s (paper -858/125) ; Li4h %s, Li5h %s, Li4h l2 %s' % (
        c['z5'], c['z2z3'], c['Li4h'], c['Li5h'], c['Li4h l2']))
    res[dps] = list(rel)
    sys.stdout.flush()
if len(res) == 2:
    a, b = list(res.values())
    print('same relation at both precisions:', a == b or a == [-x for x in b])
