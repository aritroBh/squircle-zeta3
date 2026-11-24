import sys, pickle; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as F
import phi
s1,s2,s3,s4,s5=sp.symbols('s1 s2 s3 s4 s5')
G=pickle.load(open('M35_sym.pkl','rb')); M3,M5=G[3],G[5]
s13=s4-s1-s2; s24=s5-s2-s3
P2=sp.Matrix([[-s3*s4+s1*(s3-s5), s13*s24],[s1*s3,(s1+s2)*(s2+s3)-s4*s5]])      # SS (3.19)
for a in [(1,1,1,2),(1,2,2,1),(3,1,2,5),(2,1,1,3),(2,3,2,3),(2,2,3,3)]:
    sub={s1:a[2],s2:0,s3:a[1],s4:a[0]+a[2],s5:a[1]+a[3]}
    m3=M3.subs(sub); m5=M5.subs(sub); p2=P2.subs(sub)
    comm=(m5*m3-m3*m5)[0,0]
    c_z3z5=(m5*m3)[0,0]+comm/5          # SS basis -> std basis (zeta^{SS}_{3,5} = z3z5 - z8 - zeta(3,5))
    c_z2z3sq=sp.Rational(1,2)*(p2*m3*m3)[0,0]
    c_35=-comm/5
    eng=phi.E.reduced_coeffs(a,8)[8]
    print(a,'SS-formula: z3z5=%s z2z3^2=%s zeta(3,5)=%s | engine: z3z5=%s z2z3^2=%s zeta(3,5)=%s'%(c_z3z5,c_z2z3sq,c_35,eng[1],eng[2],eng[3]),
          ' ALL MATCH:',sp.Rational(eng[1].numerator,eng[1].denominator)==c_z3z5 and sp.Rational(eng[2].numerator,eng[2].denominator)==c_z2z3sq and sp.Rational(eng[3].numerator,eng[3].denominator)==c_35)
# weight-2 check of P2_11 vs engine zeta(2) coefficient
for a in [(1,1,1,2),(3,1,2,5)]:
    sub={s1:a[2],s2:0,s3:a[1],s4:a[0]+a[2],s5:a[1]+a[3]}
    print(a,'P2_11 =',P2.subs(sub)[0,0],' engine zeta2 coeff =',phi.E.reduced_coeffs(a,2)[2][0])
