import pickle, sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
H3=sp.sympify(pickle.load(open('H7_div_C4_quotient.pkl','rb')))
print("H3 (deg-3 control quotient) =", sp.expand(H3))
print("H3 has all-positive coefficients:", all(c>0 for c in sp.Poly(sp.expand(H3),*A).coeffs()))
P7=sp.sympify(open('P7.txt').read())
P7p=sp.Poly(P7,*A)
print("\nP7: terms", len(sp.Add.make_args(P7)), "deg", P7p.total_degree(), "content", sp.gcd(list(P7p.coeffs())))
print("P7(1,1,1,1,1) =", P7.subs(dict(zip(A,[1,1,1,1,1]))))
print("P7 coeff sign pattern: positive terms =", sum(1 for c in P7p.coeffs() if c>0), "negative =", sum(1 for c in P7p.coeffs() if c<0))
