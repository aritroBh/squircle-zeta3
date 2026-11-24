import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
t1 = 1/(64*m**2*(2*m+1)**4)
t2 = 1/(64*(m+1)**2*(2*m+1)**4)
S1 = sp.summation(t1, (m, 1, sp.oo))
S2 = sp.summation(t2, (m, 1, sp.oo))
Sigma = sp.Rational(1,64) + S1 + S2
print("Sigma =", sp.simplify(Sigma))
print(sp.N(Sigma, 40))
# then b1''
N2 = sp.log(2)**2 - 3*sp.log(2) + 5*sp.pi**2/24
b1p = (sp.log(2)-1)/2
b1pp = sp.simplify(N2/2 - b1p**2 - Sigma)
print("b1'' =", b1pp)
print(sp.N(b1pp,40))
