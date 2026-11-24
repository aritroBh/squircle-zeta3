
import pickle, sympy as sp
from fractions import Fraction as F
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
forms4 = {k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open(f'{chn1}/forms.pkl','rb')).items()}
x,y,z,w = sp.symbols('x y z w')
Q5_4 = sp.cancel(sp.Rational(224)*sum(sp.Rational(c.numerator,c.denominator)*x**(e[0]-1)*y**(e[1]-1)*z**(e[2]-1)*w**(e[3]-1) for e,c in forms4[(10,6)].items())/(x+y+z+w))
Q5_4 = sp.expand(Q5_4)
print("Q5(P4) terms:", len(sp.Poly(Q5_4).terms()), "leading:", sp.Poly(Q5_4).terms()[:3])

a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV=[a1,a2,a3,a4,a5]
forms10 = pickle.load(open('forms10.pkl','rb'))
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
Sig = a1+a2+a3+a4+a5
P10_z37 = to_sym(forms10[10][6])
Q7 = sp.expand(sp.div(P10_z37, a2*a3*a4*Sig, *AV)[0])
# content of Q7
from math import gcd
from functools import reduce
coeffs=[t[1] for t in sp.Poly(Q7,*AV).terms()]
g = reduce(gcd, [abs(sp.Rational(c).p) for c in coeffs])
den = 1
for c in coeffs:
    den = sp.ilcm(den, sp.Rational(c).q)
print("Q7 coeffs: numerator gcd", g, "denominator lcm", den)
r0 = sp.expand(Q7.subs(a1,0)*224 - a5*Q5_4.subs({x:a2,y:a3,z:a4,w:a5}))
print("224*Q7(a1=0) == a5*Q5(P4)(a2..a5):", r0==0)
# also a5=0 side
r5 = sp.expand(Q7.subs(a5,0)*224 - a1*Q5_4.subs({x:a1,y:a2,z:a3,w:a4}))
print("224*Q7(a5=0) == a1*Q5(P4)(a1..a4):", r5==0)
