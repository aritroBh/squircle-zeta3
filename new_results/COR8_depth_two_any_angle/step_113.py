
import sys, pickle, re, itertools, time
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 90
from fractions import Fraction as Fr
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
# raise K
import cor8_series as CS
CS.K = 120
DS = {eval(k): mp.mpf(v) for k,v in pickle.load(open('series_D_v3.pkl','rb'))['DS'].items()}
def sLi_n(g, n):
    g0 = g.subs(c, 0)
    dlogg = sdiff(slog_of_rat(g, c))
    if g0 == 1:
        cur = sscale(slog_of_rat(1 - g, c), -1)
        for k in range(2, n+1):
            cur = sintegrate(smul(cur, dlogg))
            cur[(0, 0)] = mp.zeta(k)
        return trim(cur)
    cur = sscale(slog_of_rat(1 - g, c), -1)
    for k in range(2, n+1):
        cur = sintegrate(smul(cur, dlogg))
        cur[(0, 0)] = mp.re(mp.polylog(k, mp.mpf(g0))) if g0 != 0 else mp.mpf(0)
    return trim(cur)
# args used in the v3 identity
S = pickle.load(open('full_sols_v3.pkl','rb'))
sols, alln = S['sols'], S['alln']
li4_args = [nm for j,(k_,nm) in enumerate(alln) if k_=='Li4' and any(sols[E][j] for E in (-1,0,1))]
print(len(li4_args), "Li4 args used")
L2 = mp.log(2)
def slog_h(h):
    return {'c': {(0,1): mp.mpf(1)}, '1-c': slog_of_rat(1-c, c), '1+c': slog_of_rat(1+c, c)}[h]
t0 = time.time()
basis = []  # (name, series)
for nm in li4_args:
    g = sp.sympify(nm)
    S3 = sLi_n(g, 3); S2 = sLi_n(g, 2)
    basis.append((f'l2*Li3({nm})', sscale(S3, L2)))
    basis.append((f'z2*Li2({nm})', sscale(S2, mp.zeta(2))))
    basis.append((f'l2^2*Li2({nm})', sscale(S2, L2**2)))
    for h in ['c','1-c','1+c']:
        basis.append((f'l2*Li2({nm})*log({h})', sscale(smul(S2, slog_h(h)), L2)))
print("arg basis done", time.time()-t0)
for h1, h2 in itertools.combinations_with_replacement(['c','1-c','1+c'], 2):
    P = smul(slog_h(h1), slog_h(h2))
    basis.append((f'z2*log({h1})*log({h2})', sscale(P, mp.zeta(2))))
    basis.append((f'l2^2*log({h1})*log({h2})', sscale(P, L2**2)))
for trip in itertools.combinations_with_replacement(['c','1-c','1+c'], 3):
    P = {(0,0): mp.mpf(1)}
    for t in trip:
        P = smul(P, slog_h(t))
    basis.append(('l2*' + '*'.join(f'log({t})' for t in trip), sscale(P, L2)))
for g, h in [('c','c'), ('1-c','1-c'), ('1+c','1+c')]:
    basis.append((f'z3*log({g})', sscale(slog_h(g), mp.zeta(3))))
    basis.append((f'z2*l2*log({g})', sscale(slog_h(g), mp.zeta(2)*L2)))
    basis.append((f'l2^3*log({g})', sscale(slog_h(g), L2**3)))
for nm, val in [('z4', mp.zeta(4)), ('z3*l2', mp.zeta(3)*L2), ('z2*l2^2', mp.zeta(2)*L2**2), ('l2^4', L2**4)]:
    basis.append((nm, {(0, 0): val}))
print("total basis:", len(basis), time.time()-t0)
pickle.dump({'basis': [(nm, {str(k): str(v) for k, v in S_.items()}) for nm, S_ in basis]}, open('invisible_series.pkl','wb'))
