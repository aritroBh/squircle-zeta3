import sys, pickle
base = '<project>/research/08_theory/'
sys.path.insert(0, base+'TOW3_tower_depth6')
import sympy as sp
from collections import Counter
L = pickle.load(open(base+'TOW3_tower_depth6/L_parsed.pkl','rb'))
d, m = sp.symbols('d m')
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
def bell(delta, L):
    tot = sp.Integer(0)
    for part in parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (L[sz]**mult) / sp.factorial(mult)
        tot += coef
    return tot
a5b = sp.cancel(bell(5, L))
print('Bell a5 at (5,3):', a5b.subs({d:5,m:3}))
print('claimed a5 at (5,3): 3/56 - 31/70 i')
print('difference:', sp.expand(a5b.subs({d:5,m:3}) - (sp.Rational(3,56)-sp.Rational(31,70)*sp.I)))
# is the Bell a5 identically equal to the claimed closed form?
N5 = ( -16*sp.I*d**7 + 304*sp.I*d**6 + (-80-2520*sp.I)*d**5 + 160*d**4*m + (1040+11840*sp.I)*d**4
      - 2080*d**3*m + (-5140-34069*sp.I)*d**3 + (10280-240*sp.I)*d**2*m + (11440+60081*sp.I)*d**2
      + 240*sp.I*d*m**2 + (-22880+960*sp.I)*d*m + (-9600-59920*sp.I)*d - 960*sp.I*m**2 + 19200*m + 25200*sp.I )
a5c = N5/(480*(2*d-5)*(2*d-3))
print('Bell a5 == claimed closed form symbolically:', sp.simplify(a5b - a5c)==0)
