import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp, mpmath as mp
mp.mp.dps=50
sqrt2 = sp.sqrt(2)
t = sp.symbols('t')
# |y|^2 = (10-7r2) t^2 + (7r2-10) t + (3-2r2)
a2 = 10-7*sqrt2; a1 = 7*sqrt2-10; a0 = 3-2*sqrt2
disc = sp.expand(a1**2 - 4*a2*a0)
print("discriminant of |y|^2:", sp.nsimplify(disc), "=", sp.N(disc))
# roots
r1 = (-a1 + sp.sqrt(disc))/(2*a2); r2 = (-a1 - sp.sqrt(disc))/(2*a2)
r1s = sp.simplify(r1); 
print("root+ :", r1s)
# numerical
n1 = complex(sp.N(r1,50)); 
print("numeric:", n1)
modsq = sp.simplify(r1*sp.conjugate(r1))
print("|t*|^2 exact:", sp.nsimplify(modsq), "=", sp.N(modsq,30))
# argument
arg = mp.atan2(n1.imag, n1.real)
print("arg(t*) =", mp.nstr(arg,30), "; arg/pi =", mp.nstr(arg/mp.pi,30))
print("3pi/8 =", mp.nstr(3*mp.pi/8,30))
print("|t*| =", mp.nstr(abs(n1),30), "; |t*|^2=", mp.nstr(abs(n1)**2,30))
# is arg exactly 3pi/8? tan(arg) = imag/real
print("tan(arg) =", mp.nstr(n1.imag/n1.real,20), " 1+sqrt2 =", mp.nstr(1+mp.sqrt(2),20))
# exact real and imaginary parts
re_t = sp.simplify((r1+sp.conjugate(r1))/2)
print("Re(t*) =", sp.nsimplify(re_t))
# imag^2 = modsq - re^2
im2 = sp.simplify(modsq - re_t**2)
print("Im(t*)^2 =", sp.nsimplify(im2), "=", sp.N(im2,20))
print("Im(t*) = sqrt(", sp.nsimplify(im2), ")")
log("step: complex letter t* has arg=3pi/8 (tan=1+sqrt2), |t*|^2=1+sqrt2/2; Re=1/2")
