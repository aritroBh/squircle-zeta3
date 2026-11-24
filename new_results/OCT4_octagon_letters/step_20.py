import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
sqrt2 = sp.sqrt(2)
t, d = sp.symbols('t d')
normals = [(sqrt2/2,sqrt2/2),(0,1),(-sqrt2/2,sqrt2/2),(-1,0)]
V0 = sp.Matrix([0, sqrt2-1]); V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
y = V0 + t*(V1-V0); yp = V1-V0
s = [sp.expand(sp.Matrix(n).dot(y)) for n in normals]
a = [s[0], s[1], s[2], -s[3]]   # positive on (0,1)
ap = [sp.diff(ai,t) for ai in a]
L1,L2,L3,L4 = sp.symbols('L1 L2 L3 L4'); L=[L1,L2,L3,L4]

# F = sum a_j e^{d L_j} = 1 + d A1 + d^2 A2 + d^3 A3 + ...
A1 = sum(a[j]*L[j] for j in range(4))
A2 = sum(a[j]*L[j]**2 for j in range(4))/2
A3 = sum(a[j]*L[j]**3 for j in range(4))/6
# log F: F=1+u, u=dA1+d^2A2+d^3A3; logF = u - u^2/2 + u^3/3
# 1/F = 1 - u + u^2 - u^3 (to d^3)
# Phi = [d B1 + d^2 B2 + d^3 B3]/F  where Bk = sum a_j' L_j^k / k!
B1 = sum(ap[j]*L[j] for j in range(4))
B2 = sum(ap[j]*L[j]**2 for j in range(4))/2
B3 = sum(ap[j]*L[j]**3 for j in range(4))/6
# Phi = (dB1+d^2B2+d^3B3)*(1 - (dA1+d^2A2) + (dA1)^2 + ...) to order d^3
Ph1 = B1
Ph2 = B2 - B1*A1
Ph3 = B3 - B2*A1 + B1*(A1**2 - A2)
Ph1=sp.expand(Ph1); Ph2=sp.expand(Ph2); Ph3=sp.expand(Ph3)
yy = y.dot(y); yyp = y.dot(yp); ypyp = yp.dot(yp)
# Q = ypyp - 2 yyp Phi + yy Phi^2 ; u = (Q - ypyp)/ypyp
u1 = sp.expand(-2*yyp*Ph1/ypyp)
u2 = sp.expand((-2*yyp*Ph2 + yy*Ph1**2)/ypyp)
u3 = sp.expand((-2*yyp*Ph3 + 2*yy*Ph1*Ph2)/ypyp)
# E3 = |y'| [u3/2 - u1 u2/4 + u1^3/16]
E3 = sp.expand(ypyp**sp.Rational(1,2)*(u3/2 - u1*u2/4 + u1**3/16))
# denominator in t
num, den = sp.fraction(sp.together(E3/ sp.sqrt(ypyp)))
print("ypyp =", sp.nsimplify(ypyp), "(constant)")
print("E3/sqrt(ypyp) denominator:", sp.factor(den))
print("denominator has t?", den.has(t))
print("max power of t in numerator:", sp.Poly(num, t).degree() if num.has(t) else 0)
log("step: E3(t) via closed form; t-denominator = %s (constant => only real log singularities)"%sp.factor(den))
