
import re
p = '<project>/research/08_theory/COR8_depth_two_any_angle/cor8_series.py'
s = open(p).read()
s = s.replace("    out = {(0, 0): mp.log(g0)}", "    out = {(0, 0): mp.log(abs(g0))}")
open(p, 'w').write(s)
print("patched")
# recompute D
import sys, pickle
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 100
from fractions import Fraction as Fr
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
D0 = pickle.load(open('series_om_dt.pkl','rb'))
OMS = {k: mp.mpf(v) for k,v in D0['OMS'].items()}
DTS = {k: mp.mpf(v) for k,v in pickle.load(open('series_dt.pkl','rb'))['DTS'].items()}
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
def s_prod(parts):
    S = {(0,0): mp.mpf(1)}
    for kind, arg in parts:
        if kind == 'log':
            S = smul(S, {'c': {(0,1): mp.mpf(1)}, '1-c': slog_of_rat(1-c, c), '1+c': slog_of_rat(1+c, c)}[arg])
        else:
            S = smul(S, sLi_n(sp.sympify(arg), int(kind[2])))
    return trim(S)
S = pickle.load(open('full_sols_v3.pkl','rb'))
sols, alln = S['sols'], S['alln']
def parse_prod(nm):
    parts = []
    i = 0
    while i < len(nm):
        m2 = re.match(r'(Li[23]|log)\(', nm[i:])
        kind = m2.group(1)
        j = i + len(kind) + 1
        depth = 1; k = j
        while depth:
            if nm[k] == '(': depth += 1
            elif nm[k] == ')': depth -= 1
            k += 1
        arg = nm[j:k-1]
        parts.append((kind, arg))
        if nm[k:k+3] == '**2':
            parts.append((kind, arg)); k += 3
        if k < len(nm):
            assert nm[k] == '*'
            k += 1
        i = k
    return parts
RHS = dict(DTS)
for j, (kind, nm) in enumerate(alln):
    cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
    if not cf: continue
    Sg = sLi_n(sp.sympify(nm), 4) if kind == 'Li4' else s_prod(parse_prod(nm))
    for E, v in cf.items():
        RHS = sadd(RHS, sscale(sshift(Sg, E), mp.mpf(v.numerator)/v.denominator))
DS = sadd(OMS, sscale(RHS, -1))
imax = max(abs(mp.im(v)) for v in DS.values())
print("max |im| in D:", mp.nstr(imax, 3))
DSr = {k: mp.re(v) for k, v in DS.items()}
items = sorted(DSr.items())
print("leading D coefficients:")
for (k, l), v in items[:30]:
    print(f"c^{k} log^{l}: {mp.nstr(v, 12)}")
pickle.dump({'DS': {(str(k)): str(v) for k,v in DSr.items()}}, open('series_D_v3.pkl','wb'))
