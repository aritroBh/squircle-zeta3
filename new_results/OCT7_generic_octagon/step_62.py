import mpmath as mp, sympy as sp, pickle
from fractions import Fraction
mp.mp.dps = 900
d = pickle.load(open('basis_indep.pkl','rb'))
names = d['names']; gen = d['gen_rels']
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2, z3 = mp.pi**2/6, mp.zeta(3)
lm = {'l2':l2,'l3':l3,'l5':l5}
import re as _re
def f2m(q): return mp.mpf(q.numerator)/q.denominator
def bval(name):
    if name=='1': return mp.mpf(1)
    if name in lm: return lm[name]
    if name=='z2': return z2
    if name=='z3': return z3
    val = mp.mpf(1); rest = name
    if rest.startswith('z2'): val*=z2; rest=rest[2:]
    mm = _re.match(r'^Li([23])\((-?\d+/\d+|-?\d+)\)', rest)
    if mm:
        q = Fraction(mm.group(2)); val *= mp.polylog(int(mm.group(1)), f2m(q)); rest = rest[mm.end():]
    for ch in _re.findall(r'l[235]', rest): val *= lm[ch]
    return val
vals900 = [bval(nm) for nm in names]
# verify 4 random relations at 900 digits and show their structure
import random
random.seed(7)
idxs = sorted(random.sample(range(len(gen)), 4))
for i in idxs:
    r = gen[i]
    res = mp.fsum([mp.mpf(c)*vals900[j] for j,c in enumerate(r)])
    terms = [(names[j], c) for j,c in enumerate(r) if c != 0]
    print(f"rel {i}: residual@900 = {mp.nstr(abs(res),3)}, {len(terms)} terms, maxcoeff {max(abs(c) for _,c in terms)}")
    print("   ", " + ".join(f"{c}*{nm}" for nm,c in terms[:10]), "..." if len(terms)>10 else "")
