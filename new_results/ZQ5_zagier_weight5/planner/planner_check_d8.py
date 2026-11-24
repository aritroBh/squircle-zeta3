"""Planner check of ZQ5's D=8 representative (out/d8/rel_min.json): T = L(chi8,5)/sqrt2 = sum c_i atom_i.
Own code: Li5 by direct series, L(chi8,5) by Hurwitz zeta; ls2 = log sqrt2, lE = log(1+sqrt2)."""
import json, re, sys
from fractions import Fraction as F
import mpmath as mp
D = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
mp.mp.dps = D + 30; EPS = mp.mpf(10) ** (-(D + 25)); r2 = mp.sqrt(2)
def fld(s):
    s = s.strip('()').replace('sqrt2', 's'); v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s', '').replace('s', ''); c = c if c not in ('', '+', '-') else c + '1'
            v += mp.mpf(F(c).numerator) / F(c).denominator * r2
        else: v += mp.mpf(F(t).numerator) / F(t).denominator
    return v
def li5(x):
    t, xn, n = mp.mpf(0), mp.mpf(1), 1
    while True:
        xn *= x; t += xn / mp.mpf(n) ** 5
        if abs(xn) < EPS: return t
        n += 1
base = {'ls2': mp.log(r2), 'lE': mp.log(1 + r2), 'pi2': mp.pi ** 2, 'z5': mp.zeta(5)}
def atom(nm):
    if nm in base: return base[nm]
    m = re.fullmatch(r'Li5\((.*)\)', nm); return li5(fld(m.group(1)))
d = json.load(open('../out/d8/rel_min.json'))
tot = mp.mpf(0); n = 0
for nm, c in zip(d['names'], d['c']):
    if nm == 'T' or F(c) == 0: continue
    v = mp.mpf(1)
    out, dep, cur = [], 0, ''
    for ch in nm:
        dep += (ch == '(') - (ch == ')')
        if ch == '*' and dep == 0: out.append(cur); cur = ''
        else: cur += ch
    for f in out + [cur]: v *= atom(f)
    tot += mp.mpf(F(c).numerator) / F(c).denominator * v; n += 1
h = [mp.zeta(5, mp.mpf(k) / 8) for k in (1, 3, 5, 7)]
L = (h[0] + h[3] - h[1] - h[2]) / mp.mpf(8) ** 5
assert F(d['c'][d['names'].index('T')]) == -1
print(f'terms {n}; RHS - L(chi8,5)/sqrt2 = {mp.nstr(tot - L / r2, 6)} at dps {D}')
