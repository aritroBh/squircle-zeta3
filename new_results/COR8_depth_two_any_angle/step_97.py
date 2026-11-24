
import sys, pickle, re
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction as Fr
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
def sLi_n(g, n):
    g0 = g.subs(c, 0)
    logg = slog_of_rat(g, c)
    dlogg = sdiff(logg)
    if g0 == 0:
        cur = sscale(slog_of_rat(1 - g, c), -1)
        for k in range(2, n+1):
            cur = sintegrate(smul(cur, dlogg))
        return trim(cur)
    if g0 == 1:
        cur = sscale(slog_of_rat(1 - g, c), -1)
        for k in range(2, n+1):
            cur = sintegrate(smul(cur, dlogg))
            cur[(0, 0)] = mp.zeta(k)
        return trim(cur)
    cur = {(0, 0): mp.re(mp.polylog(1, mp.mpf(g0)))}
    for k in range(2, n+1):
        cur = sintegrate(smul(cur, dlogg))
        cur[(0, 0)] = mp.re(mp.polylog(k, mp.mpf(g0)))
    return trim(cur)
S = pickle.load(open('full_sols_v3.pkl','rb'))
sols, alln = S['sols'], S['alln']
for j, (kind, nm) in enumerate(alln):
    cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
    if not cf or kind != 'Li4': continue
    Sg = sLi_n(sp.sympify(nm), 4)
    imax = max(abs(mp.im(v)) for v in Sg.values())
    if imax > mp.mpf('1e-40'):
        print("COMPLEX:", nm, "max im =", mp.nstr(imax, 3))
