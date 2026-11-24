import sys, random, itertools
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
from fractions import Fraction as F
import sympy as sp, engine
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
# monomials of degree 5 in 4 vars
mons=[e for e in itertools.product(range(6),repeat=4) if sum(e)==5]
print("monomials:",len(mons))
random.seed(7)
pts=[tuple(random.randint(-6,9) for _ in range(4)) for _ in range(90)]
rows=[];rhs=[]
for p in pts:
    r=engine.reduced_coeffs(tuple(F(x) for x in p),5)
    rows.append([sp.Integer(1)*sp.prod([sp.Integer(p[i])**m[i] for i in range(4)]) for m in mons])
    c=r[5][0]; rhs.append(sp.Rational(c.numerator,c.denominator))
A=sp.Matrix(rows); b=sp.Matrix(rhs)
# solve least squares exactly via normal equations (full column rank expected)
sol=(A.T*A).LUsolve(A.T*b)
res=(A*sol-b)
print("max residual:",max(abs(x) for x in res))
M5_11=sum(sol[i]*a1**mons[i][0]*a2**mons[i][1]*a3**mons[i][2]*a4**mons[i][3] for i in range(len(mons)))
M5_11=sp.expand(M5_11)
print("M5_11 on slice (zeta5 coeff of F^(23) at s23=0) =",sp.factor(M5_11))
s1,s2,s3,s4,s5=a2,0,a3,a2+a4,a1+a3
def m3_11(s1,s2,s3,s4,s5): return s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
M3_11=m3_11(s1,s2,s3,s4,s5)
M3_12=sp.expand(M3_11-s4*s5*(s4+s5))
M5_12=sp.expand(M5_11-((s4+s5)**5-s4**5-s5**5)/5)
perm={a1:a3,a2:a4,a3:a1,a4:a2}
M3_21=sp.expand(M3_12.subs(perm,simultaneous=True)); M5_21=sp.expand(M5_12.subs(perm,simultaneous=True))
# check M3_12 against SS m12 and M3_21 against m21 symbolically
m12=-(s4-s1-s2)*(s5-s2-s3)*(s1+s2+s3+s4+s5); m21=s1*s3*(s1+s2+s3-2*(s4+s5))
print("M3_12==SS m12:",sp.expand(M3_12-m12)==0," M3_21==SS m21:",sp.expand(M3_21-m21)==0)
lit=sp.expand(-(M5_12*M3_21-M3_12*M5_21)/5)
C3=lambda x1,x2,x3,x4: 2*((x1+x4)-(x2+x3))*(2*(x1+x4)+(x2+x3))*((x1+x4)+2*(x2+x3))-9*x1*x4*((x1+x4)+2*(x2+x3))+9*x2*x3*(2*(x1+x4)+(x2+x3))
P8=sp.expand(a1*a2*a3*a4*(a1+a2+a3+a4)*C3(a1,a2,a3,a4)/10)
print("POLYNOMIAL IDENTITY  -(1/5)[M5,M3]_11 == P8 :",sp.expand(lit-P8)==0)
print("factorisation of literature-side polynomial:",sp.factor(lit))
