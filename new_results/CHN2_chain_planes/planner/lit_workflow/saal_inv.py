import sys
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
import sympy as sp
from fractions import Fraction as F
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
def C3(x1,x2,x3,x4):
    p=x1+x4; q=x2+x3
    return 2*(p-q)*(2*p+q)*(p+2*q)-9*x1*x4*(p+2*q)+9*x2*x3*(2*p+q)
def P8(x1,x2,x3,x4): return x1*x2*x3*x4*(x1+x2+x3+x4)*C3(x1,x2,x3,x4)/10
ap=(a1,-a1-a3,-a1-a2,a1+a2+a3+a4)           # Saalschutz involution (nonterminating, II.24 type)
print("C3(a')-C3(a) =",sp.expand(C3(*ap)-C3(a1,a2,a3,a4)))
B0=-a2*a3/((a1+a2)*(a1+a3))
print("P8(a) + B0*P8(a') =",sp.simplify(P8(a1,a2,a3,a4)+B0*P8(*ap)))
# other choice of e (swap roles of bottom parameters): u1<->u4 reversal conjugate. check the reversed involution too
ap2=(a1+a2+a3+a4,-a4-a3,-a4-a2,a4)
print("C3(a'')-C3(a) =",sp.expand(C3(*ap2)-C3(a1,a2,a3,a4)))
# engine check at positive point -> point with negative entries
import engine
for a in [(1,1,1,2),(1,2,2,1),(3,1,2,5)]:
    A=tuple(F(x) for x in a)
    Ap=(A[0],-A[0]-A[2],-A[0]-A[1],A[0]+A[1]+A[2]+A[3])
    e=engine.reduced_coeffs(Ap,8)
    P8p_engine=e[8][3]
    P8p_formula=P8(*[sp.Rational(x.numerator,x.denominator) for x in Ap])
    P8_formula=P8(*a)
    B0v=sp.Rational(-a[1]*a[2],(a[0]+a[1])*(a[0]+a[2]))
    print("a=",a,"a'=",tuple(map(str,Ap)),"| engine zeta(3,5) coeff at a' =",P8p_engine,"| closed form at a' =",P8p_formula,"| -B0*P8(a') =",-B0v*P8p_formula,"| P8(a) =",P8_formula)
