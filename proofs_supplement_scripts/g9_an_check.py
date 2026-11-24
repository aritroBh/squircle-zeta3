# Part 3: is P^(5)(1) of the diamond = 120*sqrt2*a_5 (Proposition prop:an) exactly?
import sympy as sp
z2,z3,z5=sp.symbols('z2 z3 z5')           # zeta(2), zeta(3), zeta(5) as free symbols
pi2=6*z2                                   # the ONLY relation used: pi^2 = 6 zeta(2)
pi4=pi2**2                                 # pi^4 = 36 zeta(2)^2
# a_5 as printed in Proposition prop:an
a5=-sp.Rational(7,24)*pi2*z3+6*z3-sp.Rational(3,4)*z5-sp.Rational(2,3)*pi2+sp.Rational(11,120)*pi4
r2=sp.sqrt(2)
paper=r2*(-90*z5-210*z2*z3+396*z2**2+720*z3-480*z2)   # sentence after prop:order5
print('paper - 120*sqrt2*a5 =',sp.simplify(sp.expand(paper-120*r2*a5)))
# coefficient check omega_5 and V_5 from the proposition
t=sp.symbols('t')
V=lambda n: sp.simplify(sp.diff(sp.sqrt((1+sp.exp(4*t))/2),t,n).subs(t,0))
print('V3,V5,V7 =',V(3),V(5),V(7))
om=lambda n: sp.Rational(4*(-1)**n,n)*(1-sp.Rational(1,2**(n-1))*V(n))
print('omega3,5,7 =',om(3),om(5),om(7),' (zeta5 coeff of a5 is -3/4:',sp.Rational(-3,4)==om(5),')')
# numerical cross-check of a_5 against the definition P(1+d)/sqrt2 = sum a_n d^n (diamond perimeter from arc length)
import mpmath as mp
mp.mp.dps=60
def P(p):
    # first quadrant of |x|^p+|y|^p=1: x=cos^(2/p), y=sin^(2/p); length = int speed; P=4*quarter
    f=lambda th: (2/p)*mp.sqrt(mp.cos(th)**(4/p-2)*mp.sin(th)**2+mp.sin(th)**(4/p-2)*mp.cos(th)**2)
    return 4*mp.quad(f,[0,mp.pi/4,mp.pi/2])
h=mp.mpf(1)/400
pts=[(k*h,P(1+k*h)/mp.sqrt(2)) for k in range(-8,9)]
# least-squares-free exact polynomial interpolation of degree 16 through 17 points, read d^5 coefficient
A=mp.matrix([[x**j for j in range(17)] for x,_ in pts]); b=mp.matrix([y for _,y in pts])
c=mp.lu_solve(A,b)
z2n,z3n,z5n=mp.zeta(2),mp.zeta(3),mp.zeta(5)
a5n=-7*mp.pi**2*z3n/24+6*z3n-3*z5n/4-2*mp.pi**2/3+11*mp.pi**4/120
print('a_5 numeric from P(p) interpolation :',mp.nstr(c[5],12))
print('a_5 from the closed form            :',mp.nstr(a5n,12))
print('a_2,a_3 check:',mp.nstr(c[2],10),mp.nstr(mp.pi**2/6,10),mp.nstr(c[3],10),mp.nstr(z3n-mp.pi**2/3,10))
