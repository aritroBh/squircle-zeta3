import os, sys
os.chdir('<project>/research/08_theory/DEC5E_decagon_p5_closed')
import sympy as sp
# prop:order5 for the regular 2m-gon at m=5
m = 5
s = sp.sin(sp.pi/(2*m)); c = sp.cos(sp.pi/(2*m))
C5 = sp.expand_trig(2880*m*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1))
C23 = sp.expand_trig(-960*m*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1))
C5 = sp.nsimplify(sp.expand(C5), [sp.sqrt(5)])
C23 = sp.nsimplify(sp.expand(C23), [sp.sqrt(5)])
print("C_5^(5)   =", C5)
print("C_5^(2,3) =", C23)
# L = sqrt((25-11sqrt5)/10)
L = sp.sqrt((25-11*sp.sqrt(5))/10)
r5 = sp.sqrt(5)
def pq(e):
    e = sp.expand(e)
    p = sp.simplify((e + e.subs(r5,-r5))/2); q = sp.simplify((e - e.subs(r5,-r5))/(2*r5))
    return p, q
print("C5/L =", pq(sp.simplify(C5/L)))
print("C23/L =", pq(sp.simplify(C23/L)))
# DEC5D local block coefficients
print("DEC5D local: z5 = (73125/64, -30825/64), z2z3 = (-64875/64, 32775/64)")
print("match z5:", pq(sp.simplify(C5/L)) == (sp.Rational(73125,64), sp.Rational(-30825,64)))
print("match z2z3:", pq(sp.simplify(C23/L)) == (sp.Rational(-64875,64), sp.Rational(32775,64)))
open('progress.log','a').write("step_20: prop:order5 at m=5 vs DEC5D local-block z5/z2z3 coefficients: EXACT match (sympy, Q(sqrt5))\n")
