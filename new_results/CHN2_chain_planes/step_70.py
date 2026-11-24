
import pickle
from fractions import Fraction as F
import sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV = [a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)

P8 = to_sym(forms10[8][3])
for i,v in enumerate(AV):
    print(f"  a{i+1} divides P8:", sp.expand(P8.subs(v,0))==0)
Sig = a1+a2+a3+a4+a5
print("Sigma divides P8:", sp.expand(sp.rem(P8, Sig, *AV))==0)
Q5 = sp.expand(P8/(a2*a3*a4))
print("Q5 poly terms:", len(sp.Poly(Q5,*AV).terms()))
print("Sigma divides Q5:", sp.expand(sp.rem(Q5, Sig, *AV))==0)
Q4p = sp.expand(Q5/Sig)
print("Q5/Sig poly terms:", len(sp.Poly(Q4p,*AV).terms()))
