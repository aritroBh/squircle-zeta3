import pickle, sympy as sp
forms=pickle.load(open('forms10.pkl','rb'))
import mzv_exact as M
f10=forms[10]
print("forms[10] type:", type(f10), "len:", len(f10))
print("W10 basis:")
for i,mon in enumerate(M.BASIS[10]):
    print(f"  [{i}] {'*'.join(','.join(map(str,w)) for w in mon)}")
# the zeta(3,7) coordinate is index 6
z37 = f10[6]
print("\nforms[10][6] (zeta(3,7) coord): type", type(z37))
if isinstance(z37, dict):
    print("terms:", len(z37), "sample:", list(z37.items())[:3])
    # convert to sympy
    a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
    z37sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in z37.items())
    Q7=sp.sympify(open('Q7p.txt').read())
    # check divisibility by Q7'
    q,r=sp.div(z37sp, Q7, a1,a2,a3,a4,a5)
    print("\nzeta(3,7) coord / Q7': remainder == 0?", sp.expand(r)==0)
    if sp.expand(r)==0:
        print("quotient:", sp.factor(q))
