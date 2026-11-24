"""Planner check of ZQ5's D=12 relation (out/d12/rel_min_atom.json): T = L(chi12,5)/sqrt3 = sum c_i atom_i.
Own code: Li5 by direct series (mp.polylog only for |x| > 0.9), L(chi12,5) by Hurwitz zeta.
Letters: lE3 = log(2+sqrt3), ls3 = log sqrt3, lu = log(1+sqrt3)."""
import json, re, sys
from fractions import Fraction as F
import mpmath as mp
D = int(sys.argv[1]) if len(sys.argv) > 1 else 600
mp.mp.dps = D + 30; EPS = mp.mpf(10) ** (-(D + 25)); r3 = mp.sqrt(3)
def fld(s):
    s = s.strip('()').replace('sqrt3', 's'); v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s', '').replace('s', ''); c = c if c not in ('', '+', '-') else c + '1'
            v += mp.mpf(F(c).numerator) / F(c).denominator * r3
        else: v += mp.mpf(F(t).numerator) / F(t).denominator
    return v
def li5(x):
    if abs(x) > 0.9: return mp.polylog(5, x)
    t, xn, n = mp.mpf(0), mp.mpf(1), 1
    while True:
        xn *= x; t += xn / mp.mpf(n) ** 5
        if abs(xn) < EPS: return t
        n += 1
base = {'lE3': mp.log(2 + r3), 'ls3': mp.log(r3), 'lu': mp.log(1 + r3), 'pi2': mp.pi ** 2, 'z5': mp.zeta(5)}
def split(nm):
    out, dep, cur = [], 0, ''
    for ch in nm:
        dep += (ch == '(') - (ch == ')')
        if ch == '*' and dep == 0: out.append(cur); cur = ''
        else: cur += ch
    return out + [cur]
def atom(nm):
    if nm in base: return base[nm]
    return li5(fld(re.fullmatch(r'Li5\((.*)\)', nm).group(1)))
d = json.load(open('../out/d12/rel_min_atom.json'))
tot = mp.mpf(0); n = 0; nli = 0
for nm, c in zip(d['names'], d['c']):
    if nm == 'T' or F(c) == 0: continue
    v = mp.mpf(1)
    for f in split(nm): v *= atom(f)
    tot += mp.mpf(F(c).numerator) / F(c).denominator * v; n += 1; nli += nm.startswith('Li5')
chi = {1: 1, 5: -1, 7: -1, 11: 1}
L = sum(chi[a] * mp.zeta(5, mp.mpf(a) / 12) for a in chi) / mp.mpf(12) ** 5
assert F(d['c'][d['names'].index('T')]) == -1
print(f'terms {n} (Li5 {nli}); L(chi12,5) = {mp.nstr(L, 25)}; RHS - L/sqrt3 = {mp.nstr(tot - L / r3, 6)} at dps {D}')
print('dps in force:', mp.mp.dps, ' RHS =', mp.nstr(tot, 40), ' L/sqrt3 =', mp.nstr(L / r3, 40))
print('raw difference (mpf):', repr(tot - L / r3))
bad = tot + mp.mpf(1) / 10**40 * atom('z5')
print('perturbation control (should be ~1e-40):', mp.nstr(bad - L / r3, 6))
