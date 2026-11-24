
import sys
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 60
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
def sLi_n(g, n):
    g0 = g.subs(c, 0)
    dlogg = sdiff(slog_of_rat(g, c))
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
def evalS(S, cc):
    lc = mp.log(cc)
    return sum(v * cc**k * lc**l for (k,l), v in S.items())
cc = mp.mpf('3e-3')
bad = 0
for nm in ['c', '-c', '1 - c', '1/(1 - c)', 'c + 1', '2/(1 - c)', '2/(c + 1)', '(c - 1)/(c + 1)',
           '1/2 - c/2', 'c/2 + 1/2', '(1 - c)*(c + 1)', '(c + 1)**2/(1 - c)**2', '1/((1 - c)*(c + 1))',
           '-c/(1 - c)', '-2*c/(1 - c)', '2*c/(c + 1)', '4*c/(c + 1)**2', '-4*c/(1 - c)**2', 'c**2',
           'c/(c + 1)', '(1 - c)/(c + 1)', '(1 - c)**2/(c + 1)**2', '-c**2/((1 - c)*(c + 1))']:
    g = sp.sympify(nm)
    d = abs(evalS(sLi_n(g, 4), cc) - mp.re(mp.polylog(4, mp.mpf(sp.N(g.subs(c, cc), 40)))))
    if d > mp.mpf('1e-40'):
        print(f"BAD {nm}: |d|={mp.nstr(d,3)}"); bad += 1
print("bad:", bad)
