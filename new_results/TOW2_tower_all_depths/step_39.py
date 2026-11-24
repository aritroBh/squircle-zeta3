import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
import sympy as sp
from fractions import Fraction as Fr

d, m = sp.symbols('d m')
i = sp.I
a2 = (-4*d**3 + 24*d**2 + d*(-49 + 2*i) - 4*i*m + 33)/(4*(2*d - 3))
a3 = (4*i*d**4 - 36*i*d**3 + d**2*(6 + 125*i) - 12*d*m + d*(-18 - 198*i) + 36*m + 117*i)/(24*d - 36)
a4 = (16*d**6 - 240*d**5 + d**4*(1528 - 48*i) + 96*i*d**3*m + d**3*(-5296 + 432*i) - 864*i*d**2*m + d**2*(10525 - 1308*i) + d*m*(48 + 2616*i) + d*(-11321 + 1320*i) - 48*m**2 - 2640*i*m + 5040)/(384*d**2 - 1536*d + 1440)
# TOW1 R2 for cross-check
def rfun(dd, mm):
    return 2*(dd-2*mm) - 2*dd*i - 4*(dd-1)*(dd-2*mm)/(2*dd-3)
R2 = sp.Rational(1,2)*((3*d-3)/2 - i*rfun(d,m)/2 - (d-2)**2)
print('a2 == R2:', sp.simplify(a2 - R2) == 0)

def aQSI(expr, dd, mm):
    val = expr.subs({d: dd, m: mm})
    val = sp.expand(val)
    # convert to QSI: val in Q(i)
    re = sp.re(val); im = sp.im(val)
    return QSI(Fr(sp.Rational(re)), 0, Fr(sp.Rational(im)), 0)

bad = 0; tested = 0
for dd in range(2, 12):
    for key in keys_of_degree(dd):
        p, q = key
        rho0 = rho(key, 0)
        for delta, af in ((3, a3), (4, a4)):
            if delta > 2*dd-4: continue
            tested += 1
            expect = rho0 * aQSI(af, dd, p)
            actual = rho(key, delta)
            if expect != actual:
                bad += 1
                if bad < 5: print('MISMATCH', key, delta, expect, actual)
print(f'rho(3), rho(4) closed forms vs initial data: {tested} cells, {bad} bad')
# also verify rho(2) form matches rho_c (TOW1) everywhere
bad2 = 0
for dd in range(2, 12):
    for key in keys_of_degree(dd):
        p,q = key
        if rho(key,2) != rho(key,0)*aQSI(a2,dd,p): bad2 += 1
print('rho(2) recheck:', bad2, 'bad')
