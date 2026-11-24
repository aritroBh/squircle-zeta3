
import pickle, sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV=[a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
P8 = to_sym(forms10[8][3])
Sig = a1+a2+a3+a4+a5
quo, rem = sp.div(P8, a2*a3*a4*Sig, *AV)
assert sp.expand(rem)==0
Q4p = sp.expand(quo)
print("Q4' degree-4 terms:", len(sp.Poly(Q4p,*AV).terms()))
# CHN1 C3 in (x,y,z,w)
x,y,z,w = sp.symbols('x y z w')
C3 = (4*(x**3+w**3) - 4*(y**3+z**3) + 3*w*x*(x+w) - 3*y*z*(y+z)
      + 6*(x**2+w**2)*(y+z) - 6*(y**2+z**2)*(x+w) + 6*y*z*(x+w) - 6*w*x*(y+z))
# restriction a1=0: Q4'(0,a2,a3,a4,a5) == a5*C3(a2,a3,a4,a5)/10 ?
r0 = sp.expand(Q4p.subs(a1,0) - a5*C3.subs({x:a2,y:a3,z:a4,w:a5})/10)
print("Q4'(a1=0) == a5*C3(a2..a5)/10:", r0==0)
r5 = sp.expand(Q4p.subs(a5,0) - a1*C3.subs({x:a1,y:a2,z:a3,w:a4})/10)
print("Q4'(a5=0) == a1*C3(a1..a4)/10:", r5==0)
# factor Q4'
f = sp.factor(Q4p)
print("factor(Q4'):", f)
