# Main-session planner check of M4c: evaluate R5 from closed_form.json with an independent parser/evaluator (mpmath only)
# and compare with the stored high-precision W5. Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m / (n^s m^t).
import json, re
from fractions import Fraction
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr
mp.dps = 80
d = json.load(open('closed_form.json'))
s2 = sqrt(2); b = s2 - 1
def ev(expr):
    return mp.mpf(eval(expr.replace('sqrt2', 's2'), {'s2': s2, '__builtins__': {}}))
def li_st(s, t, x, y, N=None):
    tot = mpf(0); inner = mpf(0); xn = mpf(1); n = 1
    while True:
        xn *= x
        term = xn * inner / mpf(n)**s
        tot += term
        if n > 5 and abs(xn) < mpf(10)**(-mp.dps - 5): break
        inner += y**n / mpf(n)**t
        n += 1
    return tot
def split_top(k):
    out, depth, cur = [], 0, ''
    for ch in k:
        if ch == '(': depth += 1
        if ch == ')': depth -= 1
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out
def factor(f):
    base = {'lam': log(2), 'L': log(1 + s2), 'pi2': pi**2, 'z3': zeta(3), 'z5': zeta(5), '1': mpf(1)}
    if f in base: return base[f]
    if f == 'Li2(b)': return polylog(2, b)
    m = re.fullmatch(r'Li32\((.*),(-?1)\)', f)
    if m: return li_st(3, 2, ev(m.group(1)), mpf(int(m.group(2))))
    m = re.fullmatch(r'Li(\d)\((.*)\)', f)
    if m: return polylog(int(m.group(1)), ev(m.group(2)))
    raise ValueError(f)
R5 = mpf(0)
for key, (cr, cs) in d['terms'].items():
    val = mpf(1)
    for f in split_top(key): val *= factor(f)
    c = mpf(Fraction(cr).numerator) / Fraction(cr).denominator + s2 * mpf(Fraction(cs).numerator) / Fraction(cs).denominator
    R5 += c * val
W5 = mpf(open('W5_4200.txt').read().split()[0])
print("R5 =", nstr(R5, 40)); print("W5 =", nstr(W5, 40)); print("R5 - W5 =", nstr(R5 - W5, 5))
