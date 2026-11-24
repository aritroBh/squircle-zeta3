from fractions import Fraction as Fr
from math import gcd
# Fast FT: element = tuple (A,B,C,d), value (A + B*al + C*al^2)/d, al^3 = -al^2+2al+1
def ft(A=0,B=0,C=0,d=1): return (Fr(A),Fr(B),Fr(C),int(d))
ZERO=(Fr(0),)*3+(1,); ONE=(Fr(1),Fr(0),Fr(0),1); AL=(Fr(0),Fr(1),Fr(0),1)
def fadd(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y
    from math import lcm
    d=lcm(d1,d2)
    return (A1*(d//d1)+A2*(d//d2), B1*(d//d1)+B2*(d//d2), C1*(d//d1)+C2*(d//d2), d)
def fneg(x): A,B,C,d=x; return (-A,-B,-C,d)
def fsub(x,y): return fadd(x,fneg(y))
def fmul(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y
    # al^3=-al^2+2al+1 ; al^4=3al^2+al-1
    A=A1*A2+(B1*C2+C1*B2)-C1*C2
    B=A1*B2+B1*A2+2*(B1*C2+C1*B2)-C1*C2
    C=A1*C2+C1*A2+B1*B2-(B1*C2+C1*B2)+3*C1*C2
    return (A,B,C,d1*d2)
def fs(x,q):  # scale by rational
    A,B,C,d=x; q=Fr(q); return (A*q.numerator,B*q.numerator,C*q.numerator,d*q.denominator)
def finv(x):
    A,B,C,d=x
    # mult matrix of (A+B al+C al^2): cols = images of 1,al,al^2 (as (const,al,al^2) coords)
    # al*(A+B al+C al^2) = C + (A+2C)al + (B-C)al^2
    # al^2*(..) = al*(above): al*(C+(A+2C)al+(B-C)al^2) = (B-C) + (C+2(B-C))al + (A+2C-(B-C))al^2 = (B-C)+(2B-C)al+(A-B+3C)al^2
    M=[[A,   C,     B-C  ],
       [B,   A+2*C, 2*B-C],
       [C,   B-C,   A-B+3*C]]
    # solve M * (a,b,c)^T = (d,0,0)^T  -> inverse image of d*1
    import sympy as sp
    Ms=sp.Matrix([[sp.Rational(v) for v in row] for row in M])
    col=Ms.inv()*sp.Matrix([sp.Rational(d),0,0])
    return (Fr(col[0]),Fr(col[1]),Fr(col[2]),1)
def fdiv(x,y): return fmul(x,finv(y))
def feq0(x): return x[0]==0 and x[1]==0 and x[2]==0
import mpmath as mp
_ALR=2*mp.cos(2*mp.pi/7)
def fev(x,prec=60):
    old=mp.mp.dps; mp.mp.dps=max(prec,30)
    A,B,C,d=x; a=2*mp.cos(2*mp.pi/7)
    v=(mp.mpf(int(A.numerator))/A.denominator+mp.mpf(int(B.numerator))/B.denominator*a+mp.mpf(int(C.numerator))/C.denominator*a*a)/d
    mp.mp.dps=old; return v
# test
x1=ft(1,-1,0); one=ONE
print("finv test:", feq0(fsub(fmul(x1,finv(x1)),ONE)))  # x1 * x1^-1 = 1
print("x1 eval:", fev(x1,30), " (should be 1-alpha=-0.24698)")
import pickle
open('ftfast.py','w').write('')  # placeholder, will write real module next
print("ftfast OK")
