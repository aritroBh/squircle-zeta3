"""Planner check: evaluate Codex's R5 (closed_form.json) with independent mpmath code and compare it with
M4a's W5_600.txt (Kimi; a different quadrature code).

Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m / (n^s m^t), by direct double series. Li_k uses mpmath.polylog.
No imports from the M4c scripts.
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent.parent
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 640
mp.mp.dps = DPS + 20
r2 = mp.sqrt(2)
lam, L = mp.log(2), mp.log(1 + r2)

def field(s):  # '(-16+12*sqrt2)', '1/2*sqrt2', '-1', '(1/2-3/8*sqrt2)', 'b'
    s = s.strip()
    if s == 'b': return r2 - 1
    s = s.strip('()')
    s = s.replace('sqrt2', 's')
    toks = re.findall(r'[+-]?[^+-]+', s)
    v = mp.mpf(0)
    for t in toks:
        if 's' in t:
            c = t.replace('*s', '').replace('s', '')
            c = c if c not in ('', '+', '-') else c + '1'
            v += mp.mpf(F(c).numerator) / F(c).denominator * r2
        else:
            v += mp.mpf(F(t).numerator) / F(t).denominator
    return v

cache = {}
def li_st(s, t, x, y):
    tot = mp.mpf(0); H = mp.mpf(0); xn = mp.mpf(1); ym = mp.mpf(1)
    eps = mp.mpf(10) ** (-(DPS + 15)); n = 1
    while True:
        xn *= x
        term = xn * H / mp.mpf(n) ** s
        tot += term
        if n > 50 and abs(xn) < eps: break
        ym *= y
        H += ym / mp.mpf(n) ** t
        n += 1
    return tot

def atom(a):
    if a in cache: return cache[a]
    if a == 'lam': v = lam
    elif a == 'L': v = L
    elif a == 'pi2': v = mp.pi ** 2
    elif a == 'z3': v = mp.zeta(3)
    elif a == 'z5': v = mp.zeta(5)
    else:
        m = re.fullmatch(r'Li(\d)\((.*)\)', a)
        if m and m.group(1) in '2345' and len(m.group(1)) == 1:
            v = mp.polylog(int(m.group(1)), field(m.group(2)))
        else:
            m = re.fullmatch(r'Li(\d)(\d)\((.*),(.*)\)', a)
            assert m, a
            v = li_st(int(m.group(1)), int(m.group(2)), field(m.group(3)), field(m.group(4)))
    cache[a] = v
    return v

def split(name):
    out, depth, cur = [], 0, ''
    for ch in name:
        if ch == '(': depth += 1
        if ch == ')': depth -= 1
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out

terms = json.loads((HERE / 'closed_form.json').read_text())['terms']
R5 = mp.mpf(0)
for name, (p, q) in terms.items():
    c = mp.mpf(F(p).numerator) / F(p).denominator + r2 * mp.mpf(F(q).numerator) / F(q).denominator
    v = mp.mpf(1)
    for f in split(name): v *= atom(f)
    R5 += c * v
W5 = mp.mpf(re.sub(r'\s', '', (HERE.parent / 'M4a_decomposition' / 'W5_600.txt').read_text()).split('=')[-1])
print('terms', len(terms), 'dps', DPS)
print('R5 =', mp.nstr(R5, 40))
print('W5 (M4a W5_600.txt) =', mp.nstr(W5, 40))
print('R5 - W5_M4a =', mp.nstr(R5 - W5, 8))
s54 = json.loads((HERE / 'W5_5400.json').read_text())['value']
print('R5 - W5_M4c(5400 bits, truncated to dps) =', mp.nstr(R5 - mp.mpf(s54), 8))
