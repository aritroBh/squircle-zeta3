
import pickle, sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV=[a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
Sig = a1+a2+a3+a4+a5
P8 = to_sym(forms10[8][3])
Q4p = sp.expand(sp.div(P8, a2*a3*a4*Sig, *AV)[0])
P10_z235 = to_sym(forms10[10][5])
Q6 = sp.expand(sp.div(P10_z235, a2*a3*a4*Sig, *AV)[0])
chk = sp.expand(Q6 + (a2+a4)*(a1+a3+a5)*Q4p)
print("z2z35 == -a2a3a4 Sig (a2+a4)(a1+a3+a5) C4:', chk==0)
# z37
P10_z37 = to_sym(forms10[10][6])
q7, r7 = sp.div(P10_z37, a2*a3*a4*Sig, *AV)
print("z37/(a2a3a4 Sig) remainder zero:", sp.expand(r7)==0)
Q7 = sp.expand(q7)
print("Q7 degree:", sp.Poly(Q7,*AV).total_degree(), "terms:", len(sp.Poly(Q7,*AV).terms()))
f7 = sp.factor(Q7)
print("factor Q7:", str(f7)[:800])
