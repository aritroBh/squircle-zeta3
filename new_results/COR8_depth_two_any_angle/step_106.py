
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
g = 2/(1-c)
S = sLi_n(g, 4)
# reference: coefficients of Li4(2/(1-c)) via taylor of the differential equation
# dLi4/dc = Li3(g) * 1/(1-c); Li3(g)' = Li2(g)/(1-c); Li2(g)' = Li1(g)/(1-c); Li1(g) = -log(-(1+c)/(1-c)) real: -log((1+c)/(1-c)) = log((1-c)/(1+c))
# So Li1(2/(1-c)) = log((1-c)/(1+c))
S1ref = sscale(slog_of_rat((1-c)/(1+c), c), 1)
print("Li1 mine:", sorted(sLi_n(g,1).items())[:4])
print("Li1 ref :", sorted(S1ref.items())[:4])
# compare full Li4 with numeric taylor
def evalS(S, cc):
    lc = mp.log(cc)
    return sum(v * cc**k * lc**l for (k,l), v in S.items())
for cs in ['3e-3','1e-3','1e-4']:
    cc = mp.mpf(cs)
    print(cs, abs(evalS(S,cc) - mp.re(mp.polylog(4, mp.mpf(sp.N(g.subs(c,cc),40))))))
# print first coefficients of S vs exact via mp.taylor-like finite diff
print("S first:", [(k, mp.nstr(S.get((k,0),0), 10)) for k in range(5)])
