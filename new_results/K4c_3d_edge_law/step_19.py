
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
import sympy as sp
from sympy import sqrt, nsimplify, simplify

d1 = mpf('2.57845695565598804024194103069'); d2 = mpf('2.59012039714135033988480373635')
c1 = mp.cos(d1); c2 = mp.cos(d2)
print("cos(dih1) =", mp.nstr(c1,55))
print("cos(dih2) =", mp.nstr(c2,55))
# is cos(dih2) in Q(sqrt3)? use high-precision nsimplify via algdep
import sympy as sp
x = sp.nsimplify(sp.Float(str(mp.nstr(c2,50))), [sp.sqrt(3)])
print("cos(dih2) nsimplify with sqrt3:", x)
y = sp.nsimplify(sp.Float(str(mp.nstr(c2**2,50))), [sp.sqrt(3)])
print("cos^2(dih2) nsimplify with sqrt3:", y, "=", sp.N(y,25), " residual:", sp.N(y-sp.Float(str(mp.nstr(c2**2,50))),5))
# minimal polynomial
print("algdep cos(dih2):", sp.minpoly(sp.Float(str(mp.nstr(c2,50))) ) if False else sp.polys.numberfields.subfield.np)
# use algdep from sympy
from sympy import algdep
print("algdep(cos dih2, deg 4, sqrt3):", algdep(sp.Float(str(mp.nstr(c2,45))), 4, sp.sqrt(3)))
print("cos(dih1) exact: -(283+48*sqrt(3))/433 =", sp.N(-(283+48*sp.sqrt(3))/433, 30))
print("residual:", sp.N(-(283+48*sp.sqrt(3))/433 - sp.Float(str(mp.nstr(c1,50))), 5))
# also sin and tan for later use
print("sin(dih1):", mp.nstr(mp.sin(d1),40), " sin^2 nsimplify:", nsimplify(sp.Float(str(mp.nstr(mp.sin(d1)**2,45))), [sp.sqrt(3)]))
print("sin(dih2):", mp.nstr(mp.sin(d2),40), " sin^2 nsimplify:", nsimplify(sp.Float(str(mp.nstr(mp.sin(d2)**2,45))), [sp.sqrt(3)]))
