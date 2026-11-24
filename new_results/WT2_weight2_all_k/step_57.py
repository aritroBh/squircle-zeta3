import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from fastred import v
S2 = sp.sqrt(2); a = 1+S2; b = S2-1; I = sp.I
eps = sp.Symbol('eps')

g20   = (1+S2)/4 - S2/(2*(v**2+1))
g11_3 = -S2*(v**2+2*v-1)/(4*(v**2+1))
R10   = (1-S2/2)/v - 2/(v+a)
Rw2   = 1/(2*v) + 1/(2*(v-a)) - 1/(v+a) - 1/(2*(v+b))
DL1 = 1/(v-b) + 1/(v+a) - 1/v
DL2 = 1/(v-a) + 1/(v+b) - 1/v
def F(R):
    return sp.together(R.subs(v, 1/v)/v**2)

Gm20 = eps**2/(1-eps) * g20
Gm11 = -eps**2 * g20 + eps**3/(1-eps) * g11_3
Wcas = -2*Gm20*DL1 - Gm11*DL2 - F(sp.together(Gm11*DL1))

num = sp.together((1-eps)*Wcas)   # should be polynomial c2 eps^2 + c3 eps^3 in eps
num = sp.expand(num)
c2 = sp.together(num.coeff(eps, 2))
c3 = sp.together(num.coeff(eps, 3))
rem = sp.together(num - c2*eps**2 - c3*eps**3)
print("degree of num in eps <=3 and no other terms:", rem == 0)
print("c2(v) =", sp.factor(c2))
print("c2 == Rw2 :", sp.simplify(sp.radsimp(sp.together(c2 - Rw2)))==0)
print("c3(v) =", sp.factor(c3), " == 0 :", sp.simplify(sp.radsimp(sp.together(c3)))==0)
# also check no poles at +-i survive in c2 (i.e. c2 has only real simple poles)
den = sp.denom(sp.together(c2))
print("c2 denominator factors:", sp.factor(den, extension=S2))
