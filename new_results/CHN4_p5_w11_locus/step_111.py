import sympy as sp, pickle
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
P7=sp.sympify(open('P7.txt').read())
P7poly=sp.Poly(P7,*A)
# reversal of exponent tuple
def rev_e(e): return (e[4],e[3],e[2],e[1],e[0])
terms=P7poly.as_dict()
bad=0
for e,c in terms.items():
    re=rev_e(e)
    c2=terms.get(re)
    if c2!=c:
        bad+=1
        if bad<=5: print(f"  mismatch: {e} coeff {c}, reversed {re} coeff {c2}")
print("reversal-mismatched monomials in P7:", bad, "of", len(terms))
# Compare P7.txt to -10*Q11_poly
Q11=pickle.load(open('Q11_poly.pkl','rb'))
P7re=-10*sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in Q11.items())
print("P7.txt == -10*Q11 ?", sp.expand(P7-P7re)==0)
# check reversal-sym of the recompute
P7repoly=sp.Poly(sp.expand(P7re),*A)
t2=P7repoly.as_dict()
bad2=sum(1 for e,c in t2.items() if t2.get(rev_e(e))!=c)
print("reversal-mismatched monomials in -10*Q11_poly:", bad2)
