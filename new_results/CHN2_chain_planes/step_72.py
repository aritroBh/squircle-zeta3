
import pickle, sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV=[a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
Sig = a1+a2+a3+a4+a5

# W=10: index 5 = z2 z(3,5), index 6 = z(3,7)
P10_z235 = to_sym(forms10[10][5])
P10_z37  = to_sym(forms10[10][6])
for name, P in [('z2z35', P10_z235), ('z37', P10_z37)]:
    print(f"--- {name}: terms={len(sp.Poly(P,*AV).terms())}")
    divs = [f"a{i+1}" for i,v in enumerate(AV) if sp.expand(P.subs(v,0))==0]
    print("  divisible by:", divs, " Sigma:", sp.expand(sp.rem(P,Sig,*AV))==0)
# factor out common factors for z2z35
q, r = sp.div(P10_z235, a2*a3*a4*Sig, *AV)
print("z2z35/(a2a3a4 Sig) remainder zero:", sp.expand(r)==0)
Q6 = sp.expand(q)
print("Q6 (deg6) terms:", len(sp.Poly(Q6,*AV).terms()), " Sigma|Q6:", sp.expand(sp.rem(Q6,Sig,*AV))==0)
f6 = sp.factor(Q6)
print("factor Q6:", str(f6)[:600])
