import pickle, sympy as sp
from fractions import Fraction as F
Q11=pickle.load(open('Q11_poly.pkl','rb'))
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
# integer polynomial P7 = -10*Q11
P7=sum((-10*sp.Rational(c.numerator,c.denominator))*sp.prod(A[i]**e[i] for i in range(5)) for e,c in Q11.items())
P7=sp.Poly(P7, a1,a2,a3,a4,a5)
print("P7 is integer poly:", all(c.is_integer for c in P7.coeffs()))
# content
cont=sp.gcd(list(P7.coeffs()))
print("content of P7:", cont)
fl=sp.factor_list(P7.as_expr())
print("factor_list: scalar=", fl[0], " num factors=", len(fl[1]))
for f,e in fl[1]:
    print("  factor exp",e,"terms",len(sp.Add.make_args(sp.expand(f))), "total_deg", sp.Poly(f,a1,a2,a3,a4,a5).total_degree() if not f.is_number else 0)
print("\nP7 irreducible over Q:", len(fl[1])==1 and fl[1][0][1]==1)
# save Q11.txt (as Q11 = -P7/10) and P7.txt
open('Q11.txt','w').write(str(sp.expand(sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in Q11.items()))))
open('P7.txt','w').write(str(sp.expand(P7.as_expr())))
print("\nsaved Q11.txt and P7.txt")
print("Q11(1,1,1,1,1) =", float(P7.as_expr().subs(dict(zip(A,[1,1,1,1,1]))))/(-10))
print("P7(1,1,1,1,1) =", P7.as_expr().subs(dict(zip(A,[1,1,1,1,1]))))
