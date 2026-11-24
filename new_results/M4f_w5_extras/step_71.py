import json, re
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 880
r2, lam, L = mp.sqrt(2), mp.log(2), mp.log(1+mp.sqrt(2))
b = r2-1
def cf(p,q): return mp.mpf(F(p).numerator)/F(p).denominator + r2*mp.mpf(F(q).numerator)/F(q).denominator
WORK='<project>/research/08_theory/M4f_w5_extras'
d = json.loads(open(f'{WORK}/m4_terms.json').read())
vals = {'z3': mp.zeta(3), 'Li3(b)': mp.polylog(3,b), 'Li3(-b)': mp.polylog(3,-b),
        'Li3(1/sqrt2)': mp.polylog(3,1/r2), 'Li2(b)': mp.polylog(2,b)}
prod = lambda k: eval(k.replace('log2','lam').replace('pi^2','mp.pi**2'), {'lam':lam,'L':L,'mp':mp, **vals})
tot = mp.mpf(0)
for k,(p,q) in d['W34'].items():
    c = cf(p,q); v = prod(k); tot += c*v
    print(f'{k:18s} ({p},{q})  {mp.nstr(c*v, 25)}')
print('my W34 =', mp.nstr(tot, 30))
# M4c term values
m4c = {
 'z3': (7*r2/32 - 3/16 - mp.mpf(1)/8)*mp.zeta(3),
 'Li3(b)': (1+r2)*mp.polylog(3,b),
 'Li3(-b)': (1-r2)*mp.polylog(3,-b),
 'Li3(1/sqrt2)': -2*r2*mp.polylog(3,1/r2),
 'Li2(b)*L': 2*L*mp.polylog(2,b),
 'Li2(b)*log2': -r2*lam*mp.polylog(2,b),
 'pi^2*L': (r2/12 - mp.mpf(1)/8)*mp.pi**2*L,
 'pi^2*log2': (-(1+r2)/24 - mp.mpf(1)/24)*mp.pi**2*lam,
 'log2^3': ((1+r2)/24 - mp.mpf(1)/12)*lam**3,
 'log2^2*L': r2/4*lam**2*L,
 'log2*L^2': -(1+r2)/2*lam*L**2,
 'L^3': (5-r2)/6*L**3,
}
t2 = sum(m4c.values())
for k in m4c: print(f'{k:18s} M4c {mp.nstr(m4c[k],25)}')
print('M4c W34 =', mp.nstr(t2, 30))
