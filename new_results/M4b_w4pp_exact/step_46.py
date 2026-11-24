import os, json, pickle, re, sys
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 110
from fractions import Fraction as F
sys.path.insert(0, os.path.join(ROOT,'paper/squircle_note/scripts'))
from qsqrt2 import QS2 as Q

d = json.load(open('stage1_exact_cz.json'))
names = d['atom_names']
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]

def parse_qs2(s):
    m = re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)", s)
    if m:
        q = F(m.group(3)); q = q if m.group(2) == '+' else -q
        return float(F(m.group(1))) + mp.sqrt(2)*float(q)
    if s.endswith('*sqrt2'): return mp.sqrt(2)*float(F(s[:-6]))
    return mp.mpf(s)

AT = pickle.load(open('stage2_atoms.pkl','rb'))['atoms']
def atom_value(enc):
    v = mp.mpf(1)
    for a in enc:
        if a[0]=='log': v *= mp.log(parse_qs2(a[1]))
        elif a[0]=='pi2': v *= mp.pi**2
        elif a[0]=='z3': v *= mp.zeta(3)
        elif a[0]=='li': v *= mp.polylog(int(a[1]), parse_qs2(a[2]))
    return v

Rpp = mp.mpf(0)
for k,(name,enc,_) in enumerate(AT):
    coef = mp.mpf(cr[k].numerator)/cr[k].denominator + mp.sqrt(2)*mp.mpf(cs[k].numerator)/cs[k].denominator
    if coef: Rpp += coef*atom_value(enc)
W4pp = mp.mpf(json.load(open(os.path.join(ROOT,'research/08_theory/M4a_decomposition/quad_90.json')))['W4pp'])
print('R'' =', mp.nstr(Rpp, 100))
print('W4pp=', mp.nstr(W4pp, 100))
print('|diff| =', mp.nstr(abs(Rpp-W4pp), 3))
