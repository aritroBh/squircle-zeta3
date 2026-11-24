
import re
p = '<project>/research/08_theory/COR8_depth_two_any_angle/cor8_series.py'
s = open(p).read()
old_start = s.index("def slog_of_rat")
old_end = s.index("def sLi(")
new_fn = '''def _log1p(u):
    """log(1+u) for series u with positive k-powers only."""
    out = {}
    um = dict(u)
    for m in range(1, K+2):
        out = sadd(out, sscale(um, (-1)**(m+1)/mp.mpf(m)))
        um = smul(um, u)
    return trim(out)
def slog_of_rat(expr, c):
    """series of log(rat): rat = (nc[nmin]/dc[dmin]) * c^(nmin-dmin) * (1+u_num)/(1+u_den),
    = log|g0| + (nmin-dmin) log c + log(1+u_num) - log(1+u_den)."""
    import sympy as sp
    num, den = sp.fraction(sp.together(expr))
    num = sp.Poly(sp.expand(num), c); den = sp.Poly(sp.expand(den), c)
    nc = num.all_coeffs()[::-1]; dc = den.all_coeffs()[::-1]
    nmin = next(i for i, co in enumerate(nc) if co); dmin = next(i for i, co in enumerate(dc) if co)
    g0 = abs(mp.mpf(sp.Rational(nc[nmin])) / mp.mpf(sp.Rational(dc[dmin])))
    out = {(0, 0): mp.log(g0)}
    m = nmin - dmin
    if m:
        out[(0, 1)] = out.get((0, 1), 0) + m
    u_num = {(i - nmin, 0): mp.mpf(sp.Rational(co))/mp.mpf(sp.Rational(nc[nmin]))
             for i, co in enumerate(nc) if i > nmin}
    u_den = {(i - dmin, 0): mp.mpf(sp.Rational(co))/mp.mpf(sp.Rational(dc[dmin]))
             for i, co in enumerate(dc) if i > dmin}
    if u_num:
        out = sadd(out, _log1p(u_num))
    if u_den:
        out = sadd(out, sscale(_log1p(u_den), -1))
    return trim(out)
'''
s = s[:old_start] + new_fn + s[old_end:]
open(p, 'w').write(s)
print("rewritten slog_of_rat")
import sys
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 50
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
print("log(1-c):", sorted(slog_of_rat(1-c, c).items())[:5])
print("log(1/(1-c)):", sorted(slog_of_rat(1/(1-c), c).items())[:5])
print("log(-c/(1-c)):", sorted(slog_of_rat(-c/(1-c), c).items())[:5])
print("d/dc log(1/(1-c)):", sorted(sdiff(slog_of_rat(1/(1-c), c)).items())[:5])
