"""M4e: assemble M4 = R5_nice + R''_nice + W3^(4) + T and check vs M4_800.

R5_nice from lchi_final_1200.json (+ L(chi8,5) + N5). Independent mpmath evaluation
(reuses planner_eval_R5's atom/field/split helpers; L(chi8,5) from Lchi8_5.json).
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.argv = [sys.argv[0], '840']
sys.path.insert(0, str(HERE.parent/'M4c_w5_closed_form'/'planner'))
import planner_eval_R5 as E   # computes M4c R5 at import (prints), sets mp.mp.dps

r2, lam, L = E.r2, E.lam, E.L
b = r2 - 1
def fr(s): f=F(s); return mp.mpf(f.numerator)/f.denominator

# R5_nice
d = json.loads((HERE/'lchi_final_1200.json').read_text())
R5n = mp.mpf(0)
for nm,ps,qs in zip(d['names'], d['ca'], d['cb']):
    a,bb = F(ps), F(qs)
    if not (a or bb): continue
    c = fr(ps) + r2*fr(qs)
    v = mp.mpf(1)
    for f in E.split(nm): v *= E.atom(f)
    R5n += c*v
Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
R5n += (fr(d['aL']) + r2*fr(d['bL']))*Lchi
cert = json.loads((HERE/'cobracket.json').read_text())
for term in cert['correction']:
    def fld(pair):
        p, q = pair
        return mp.mpf(F(p).numerator)/F(p).denominator + r2*mp.mpf(F(q).numerator)/F(q).denominator
    R5n += fld(term['coefficient'])*E.li_st(3,2,fld(term['x']),fld(term['y']))
print('R5_nice =', mp.nstr(R5n, 40))
s54 = json.loads((HERE/'W5_5400.json').read_text())['value']
print('R5_nice - W5(5400) =', mp.nstr(R5n - mp.mpf(s54), 8))

# W4pp (R''_nice) from M4a
alias = {'log2':'lam','pi^2':'pi2','Li2(b)':'Li2(b)','Li3(b)':'Li3(b)',
         'Li3(-b)':'Li3((1-1*sqrt2))','Li3(1/sqrt2)':'Li3(1/2*sqrt2)'}
W4pp = mp.mpf(0)
for c, name in json.loads((HERE.parent/'M4a_decomposition'/'w4pp_closed_form.json').read_text())['closed_form_terms']:
    v = fr(c)
    for f in E.split(name):
        if f=='sqrt2': v *= r2; continue
        v *= E.atom(alias.get(f,f))
    W4pp += v
# W3^(4) = W3' - zeta(3)/8 - pi^2 lam/24 - lam^3/12  (W3' from thm:M3)
W3p = (7*r2/32 - mp.mpf(3)/16)*mp.zeta(3) + (1+r2)*mp.polylog(3,b) + (1-r2)*mp.polylog(3,-b) \
      - 2*r2*mp.polylog(3,1/r2) + (2*L - r2*lam)*mp.polylog(2,b) \
      + (r2/12 - mp.mpf(1)/8)*mp.pi**2*L - (1+r2)/24*mp.pi**2*lam + (1+r2)/24*lam**3 \
      + r2/4*lam**2*L - (1+r2)/2*lam*L**2 + (5-r2)/6*L**3
W34 = W3p - mp.zeta(3)/8 - mp.pi**2*lam/24 - lam**3/12
T = (mp.pi**2 + 6*lam**2)/24
M4 = R5n + W4pp + W34 + T
M4_ref = mp.mpf(re.sub(r'\s','',(HERE.parent/'M4a_decomposition'/'M4_800.txt').read_text()).split('=')[-1])
print('M4 (R5_nice assembly) =', mp.nstr(M4, 40))
print('M4 - M4_800 =', mp.nstr(M4 - M4_ref, 8))
