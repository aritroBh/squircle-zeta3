
import sympy as sp
import pickle
phi = sp.Symbol('phi')
h, Sh, Ch = sp.symbols('h Sh Ch', positive=True)
cf = sp.Symbol('cf')
data = pickle.load(open('monolist.pkl','rb'))
terms = [(b,p,q,sp.sympify(r)) for b,p,q,r in data]
print("The 6 terms (cos^b phi)(Sh^p)(Ch^q)(coeff):")
for b,p,q,r in terms:
    print(f"  b={b} p={p} q={q} coeff={sp.factor(r)}")
# Integrand term = r * Sh^p Ch^q * cos^b phi ; integral over [-h,h]
# int cos^-2 = tan ; int cos^-4 = tan + tan^3/3
def INT(b):
    if b==-2: return 2*sp.tan(h)
    if b==-4: return 2*(sp.tan(h)+sp.tan(h)**3/3)
I = 0
for b,p,q,r in terms:
    I += r*Sh**p*Ch**q*INT(b)
I = sp.expand(I)
# express tan(h)=Sh/Ch
I = I.subs(sp.tan(h), Sh/Ch)
I = sp.simplify(sp.expand(I))
print("\nI(h) = int_{-h}^{h} (G21+G12) dphi =")
print(sp.factor(I))
# Simplify with Ch^2 = 1-Sh^2
Is = sp.expand(I).subs(Ch**2, 1-Sh**2)
Is = sp.simplify(Is)
print("\nI(h) reduced (Ch^2=1-Sh^2):")
print(sp.factor(Is))
with open('I_of_h.txt','w') as f: f.write(sp.srepr(sp.factor(I))+"\n###\n"+sp.srepr(sp.factor(Is)))
print("saved I_of_h.txt")
