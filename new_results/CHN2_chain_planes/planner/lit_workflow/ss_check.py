import sys, pickle; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as F
import phi
s1,s2,s3,s4,s5=sp.symbols('s1 s2 s3 s4 s5')
Ms=pickle.load(open('M35_sym.pkl','rb'))
M3,M5=Ms[3],Ms[5]
s13=s4-s1-s2; s24=s5-s2-s3
m11=s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
m12=-s13*s24*(s1+s2+s3+s4+s5)
m21=s1*s3*(s1+s2+s3-2*(s4+s5))
m22=(s2+s3)*((s1+s2)*(s1+s3)-2*s1*s4)-(2*s1*s3-s4**2+2*s2*(s3+s4))*s5+s4*s5**2   # s_4^2 as printed in SS
SS3=sp.Matrix([[m11,m12],[m21,m22]])
print('M3(BSST associator) - M3(SS eq 3.19, with s_4^2 in m22) =',sp.expand(M3-SS3).tolist())
print('M5[0,0] =',sp.factor(M5[0,0]))
print('M5[0,1] =',sp.factor(M5[0,1]))
print('M5[1,0] =',sp.factor(M5[1,0]))
print('M5[1,1] =',sp.factor(M5[1,1]))
# SS-formula prediction for zeta(3,5)_std coefficient = -(1/5)[M5,M3]_{11}   (zeta^{SS}_{3,5} = zeta(5,3)_std = z3 z5 - z8 - zeta(3,5)_std)
Cm=sp.expand(M5*M3-M3*M5)
pred11=sp.factor(-Cm[0,0]/5)
print('-(1/5)[M5,M3]_11 =',pred11)
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
sub={s1:a3,s2:0,s3:a2,s4:a1+a3,s5:a2+a4}
P8=a1*a2*a3*a4*(a1+a2+a3+a4)*(2*(a1+a4-a2-a3)*(2*(a1+a4)+a2+a3)*((a1+a4)+2*(a2+a3))-9*a1*a4*((a1+a4)+2*(a2+a3))+9*a2*a3*(2*(a1+a4)+(a2+a3)))/10
print('symbolic: pred11|dictionary - P8 =',sp.expand(pred11.subs(sub,simultaneous=True)-P8))
for pt,val in (((1,1,1,2),58),((1,2,2,1),F(-1464,5)),((3,1,2,5),36696),((2,1,1,3),1134)):
    d=dict(zip((a1,a2,a3,a4),pt))
    print(pt,'SS-formula -(1/5)[M5,M3]_11 =',pred11.subs(sub,simultaneous=True).subs(d),' P8 =',sp.nsimplify(sp.Rational(val.numerator,val.denominator)) if isinstance(val,F) else val)
