import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
sqrt2 = sp.sqrt(2)
t, d = sp.symbols('t d')
normals = [(sqrt2/2,sqrt2/2),(0,1),(-sqrt2/2,sqrt2/2),(-1,0)]
V0 = sp.Matrix([0, sqrt2-1]); V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
y = V0 + t*(V1-V0); yp = V1-V0
s = [sp.expand(sp.Matrix(n).dot(y)) for n in normals]
a = [s[0], s[1], s[2], -s[3]]
ap = [sp.diff(ai,t) for ai in a]
L1,L2,L3,L4 = sp.symbols('L1 L2 L3 L4'); L=[L1,L2,L3,L4]
A1 = sum(a[j]*L[j] for j in range(4)); A2 = sum(a[j]*L[j]**2 for j in range(4))/2
A3 = sum(a[j]*L[j]**3 for j in range(4))/6
B1 = sum(ap[j]*L[j] for j in range(4)); B2 = sum(ap[j]*L[j]**2 for j in range(4))/2
B3 = sum(ap[j]*L[j]**3 for j in range(4))/6
Ph1=B1; Ph2=B2-B1*A1; Ph3=B3-B2*A1+B1*(A1**2-A2)
yy=y.dot(y); yyp=y.dot(yp); ypyp=yp.dot(yp)
u1=sp.expand(-2*yyp*Ph1/ypyp); u2=sp.expand((-2*yyp*Ph2+yy*Ph1**2)/ypyp)
u3=sp.expand((-2*yyp*Ph3+2*yy*Ph1*Ph2)/ypyp)
E3=sp.expand(u3/2 - u1*u2/4 + u1**3/16)   # dropped overall sqrt(ypyp) const
# monomial structure in L
P=sp.Poly(E3, L1,L2,L3,L4)
mons = P.terms()
print("total L-monomials in E3:", len(mons))
# which weight-3 L monomials appear (deg 3 in L)
w3 = [m for m in mons if sum(m[0])==3]
print("weight-3 L-monomials:", len(w3))
for m in w3:
    print("   ", m[0], end="")
print()
# does L1 L2 (both nonlocal) appear in weight 3? e.g. L1^2 L2, L1 L2^2, and mixed with L3,L4
cross12 = [m for m in mons if m[0][0]>0 and m[0][1]>0]
print("\nmonomials containing BOTH L1 and L2:", len(cross12))
for m in cross12: print("   ", m[0])
log("step: E3 L-monomial structure extracted; %d total, %d with both L1,L2"%(len(mons),len(cross12)))
