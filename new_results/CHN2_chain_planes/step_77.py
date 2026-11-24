
import pickle, sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
x,y,z,t = sp.symbols('x y z t')
AV=[a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
Sig = a1+a2+a3+a4+a5
C4 = sp.expand(sp.div(to_sym(forms10[8][3]), a2*a3*a4*Sig, *AV)[0])*10
Q7p = sp.expand(sp.div(to_sym(forms10[10][6]), a2*a3*a4*Sig, *AV)[0])*224
Q7p = sp.expand(Q7p)
print("C4 irreducible:", sp.factor(C4)==C4, " Q7' irreducible:", sp.factor(Q7p)==Q7p)
fams = {
 'P1-like (x,x,y,y,z)': (x,x,y,y,z),
 'P2-like (x,y,x,y,z)': (x,y,x,y,z),
 'palindromic (x,y,z,y,x)': (x,y,z,y,x),
 '(x,x,y,x,x)': (x,x,y,x,x),
 'revP1 (z,y,y,x,x)': (z,y,y,x,x),
 'revP2 (z,y,x,y,x)': (z,y,x,y,x),
 '(x,y,y,x,z)': (x,y,y,x,z),
 '(x,x,x,y,z)': (x,x,x,y,z),
 '(x,y,x,x,z)': (x,y,x,x,z),
 'diagonal (x,x,x,x,x)': (x,x,x,x,x),
}
for name, sub in fams.items():
    r1 = sp.expand(C4.subs(dict(zip(AV,sub))))
    r2 = sp.expand(Q7p.subs(dict(zip(AV,sub))))
    print(f"{name:28s} C4=0: {r1==0!s:5}  Q7'=0: {r2==0!s:5}   C4val={str(r1)[:60]}")
