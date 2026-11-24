"""M4f: M4 assembly check at ~850 digits, parsing m4_terms.json (the SAME coefficient data
that w5_closed_form.tex was generated from). Independent mpmath evaluation:
Li_{s,t} by direct double series, Li_k by mpmath.polylog, L(chi8,5) from Lchi8_5.json.
M4 = R5 + N5 + W4pp + W3^(4) + T  vs  M4a_decomposition/M4_800.txt (expect ~1e-800).
Also R5+N5 vs W5_5400 (sanity).
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 850
mp.mp.dps = DPS + 30
r2, lam, L = mp.sqrt(2), mp.log(2), mp.log(1+mp.sqrt(2))

def field(s):
    s = s.strip()
    if s == 'b': return r2 - 1
    if s == '-b': return 1 - r2
    if s == '1/sqrt2': return 1/r2
    s = s.strip('()').replace('sqrt2','s')
    v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s','').replace('s','')
            c = c if c not in ('','+','-') else c+'1'
            v += mp.mpf(F(c).numerator)/F(c).denominator*r2
        else:
            v += mp.mpf(F(t).numerator)/F(t).denominator
    return v

cache = {}
def li_st(s, t, x, y):
    tot, H, xn, ym = mp.mpf(0), mp.mpf(0), mp.mpf(1), mp.mpf(1)
    eps = mp.mpf(10)**(-(DPS+25)); n = 1
    while True:
        xn *= x; tot += xn*H/mp.mpf(n)**s
        if n > 50 and abs(xn) < eps: break
        ym *= y; H += ym/mp.mpf(n)**t; n += 1
    return tot

def atom(a):
    if a in cache: return cache[a]
    if a in ('lam','log2'): v = lam
    elif a == 'L': v = L
    elif a in ('pi2','pi^2'): v = mp.pi**2
    elif a == 'z3': v = mp.zeta(3)
    elif a == 'z5': v = mp.zeta(5)
    elif a == 'Lchi8_5':
        v = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    else:
        m = re.fullmatch(r'Li(\d)\((.*)\)', a)
        if m:
            v = mp.polylog(int(m.group(1)), field(m.group(2)))
        else:
            raise ValueError(a)
    cache[a] = v
    return v

def split(nm):
    out, depth, cur = [], 0, ''
    for ch in nm:
        if ch == '(' : depth += 1
        if ch == ')' : depth -= 1
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out

def cf(p, q):
    return mp.mpf(F(p).numerator)/F(p).denominator + r2*mp.mpf(F(q).numerator)/F(q).denominator

def prodval(key):
    key = key.replace('log2^3','log2*log2*log2').replace('log2^2','log2*log2').replace('L^3','L*L*L').replace('L^2','L*L')
    v = mp.mpf(1)
    for f in split(key): v *= atom(f)
    return v

d = json.loads((HERE/'m4_terms.json').read_text())
R5 = mp.mpf(0)
for nm,(p,q) in d['R5'].items():
    R5 += cf(p,q)*prodval(nm)
N5 = mp.mpf(0)
for t in d['N5']:
    x = cf(t['x'][0], t['x'][1]); y = cf(t['y'][0], t['y'][1])
    N5 += cf(t['c'][0], t['c'][1])*li_st(3,2,x,y)
W4 = mp.mpf(0)
for k,(p,q) in d['W4pp'].items():
    W4 += cf(p,q)*prodval(k)
W34 = mp.mpf(0)
for k,(p,q) in d['W34'].items():
    W34 += cf(p,q)*prodval(k)
T = mp.mpf(0)
for k,(p,q) in d['T'].items():
    T += cf(p,q)*prodval(k)
M4 = R5 + N5 + W4 + W34 + T
ref = mp.mpf(re.sub(r'\s','',(HERE.parent/'M4a_decomposition'/'M4_800.txt').read_text()).split('=')[-1])
w5ref = mp.mpf(json.loads((HERE/'W5_5400.json').read_text())['value'])
print(f'DPS={DPS}')
print('R5+N5 - W5_5400      =', mp.nstr(R5+N5 - w5ref, 8))
print('M4 (assembled)       =', mp.nstr(M4, 50))
print('M4 - M4_800.txt      =', mp.nstr(M4 - ref, 8))
