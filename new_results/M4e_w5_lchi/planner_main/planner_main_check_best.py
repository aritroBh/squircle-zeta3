# Main-session planner check of M4e: W5 = N5 (4 Li_{3,2}, from cobracket.json 'correction') + R5_nice_reduced.json terms
# (incl. L(chi8,5) via Hurwitz zeta), evaluated with my own parser (mpmath only), vs the certified W5_5400.json value.
# Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m / (n^s m^t)  (paper/OCT convention, validated on M4c to 2.3e-76).
import json, re
from fractions import Fraction
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr
mp.dps = 130
D = '<project>/research/08_theory/M4e_w5_lchi/'
d = json.load(open('<project>/research/08_theory/M4f_w5_extras/R5_best.json')); cb = json.load(open(D + 'cobracket.json'))
s2 = sqrt(2); b = s2 - 1
def ev(expr): return mp.mpf(eval(expr.replace('sqrt2', 's2'), {'s2': s2, '__builtins__': {}}))
def li_st(s, t, x, y):
    tot = mpf(0); inner = mpf(0); xn = mpf(1); n = 1
    while True:
        xn *= x; tot += xn*inner/mpf(n)**s
        if n > 5 and abs(xn) < mpf(10)**(-mp.dps-5): break
        inner += y**n/mpf(n)**t; n += 1
    return tot
def split_top(k):
    out, depth, cur = [], 0, ''
    for ch in k:
        if ch == '(': depth += 1
        if ch == ')': depth -= 1
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out
LCHI = sum(s*zeta(5, mpf(a)/8) for a, s in ((1, 1), (3, -1), (5, -1), (7, 1)))/mpf(8)**5
def factor(f):
    base = {'lam': log(2), 'L': log(1+s2), 'pi2': pi**2, 'z3': zeta(3), 'z5': zeta(5), '1': mpf(1), 'Lchi8_5': LCHI}
    if f in base: return base[f]
    if f == 'Li2(b)': return polylog(2, b)
    m = re.fullmatch(r'Li(\d)\((.*)\)', f)
    if m: return polylog(int(m.group(1)), ev(m.group(2)))
    raise ValueError(f)
def q(pair): return mpf(Fraction(pair[0]).numerator)/Fraction(pair[0]).denominator + s2*mpf(Fraction(pair[1]).numerator)/Fraction(pair[1]).denominator
R = mpf(0)
for key, cc in d['terms'].items():
    val = mpf(1)
    for f in split_top(key): val *= factor(f)
    R += q(cc)*val
N5 = mpf(0)
for t in cb['correction']:
    x = q(t['x']); y = q(t['y']); N5 += q(t['coefficient'])*li_st(3, 2, x, y)
W5 = mpf(json.load(open(D + 'W5_5400.json'))['value'][:140])
print('terms:', len(d['terms']), ' L(chi8,5) =', nstr(LCHI, 20), ' coefficient:', d['terms'].get('Lchi8_5'))
print('N5 + R =', nstr(N5 + R, 50)); print('W5     =', nstr(W5, 50)); print('diff   =', nstr(N5 + R - W5, 5))
