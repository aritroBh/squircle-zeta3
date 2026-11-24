import sys, os
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
ctx.prec = 600
Mb = arb('0','2.81e-180')
print('Mb**2:', (Mb**2).str(6), (Mb**2).is_finite())
La = arb('0','3e-180')
print('La**5:', (La**5).str(6))
x = arb(2); print('x**3', (x**3).str(6))
# maybe the issue: mag() returns acb-flavored? test types
r = acb(2).sqrt(); a = 1+r
va = a+acb(arb(0, arb(2)**-700))
def mag(z): return abs(z.real)+abs(z.imag)
La2 = mag(((va-0.414)*(va+a)/(4*va)).log())
print(type(La2), La2.str(6))
print((La2**2).str(6))
