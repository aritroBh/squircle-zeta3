import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
sqrt2 = sp.sqrt(2)
t, d = sp.symbols('t d')
# complex letter simplified
tstar = (1 + sp.I*(1+sqrt2))/2
print("t* =", tstar, " conj =", sp.conjugate(tstar))
print("Im(t*)=(1+sqrt2)/2 check:", sp.simplify(sp.im(tstar) - (1+sqrt2)/2))
print("|t*|^2 =", sp.simplify(tstar*sp.conjugate(tstar)))

# Build E3(t) symbolically: four terms, a_j linear, Lj = log a_j symbolic
normals = [(sqrt2/2,sqrt2/2),(0,1),(-sqrt2/2,sqrt2/2),(-1,0)]
V0 = sp.Matrix([0, sqrt2-1]); V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
y = V0 + t*(V1-V0); yp = V1-V0
# unsigned term values on edge: a3,a4 local (linear), a1,a2 nonlocal
# signed s_j; take absolute per sign on (0,1): s1>0,s2>0,s3>0,s4<0
s = [sp.expand(sp.Matrix(n).dot(y)) for n in normals]
a = [s[0], s[1], s[2], -s[3]]   # all positive on (0,1)
L1,L2,L3,L4 = sp.symbols('L1 L2 L3 L4')
Ls = [L1,L2,L3,L4]
F = sum(a[j]*sp.exp(d*Ls[j]) for j in range(4))
Fp = (1+d)*sum(sp.diff(a[j],t)*sp.exp(d*Ls[j]) for j in range(4))
print("F(d=0) =", sp.simplify(F.subs(d,0)))  # should be 1
Phi = Fp/((1+d)*F)
G = F**(-sp.Rational(2)/(1+d))*(yp.dot(yp) - 2*y.dot(yp)*Phi + y.dot(y)*Phi**2)
ser = sp.series(sp.sqrt(G), d, 0, 4).removeO()
E3 = sp.expand(ser.coeff(d,3))
# denominator in t
num, den = sp.fraction(sp.together(E3))
print("\nDenominator of E3(t):", sp.factor(den))
print("Is denominator t-dependent?", den.has(t))
log("step: E3(t) built; denominator = " + str(sp.factor(den)))
