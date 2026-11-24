import sympy as sp
x = sp.symbols('x')
y1, y2, y3, y4 = sp.symbols('y1 y2 y3 y4')
e = sp.log(x)*sp.log(x+1)/(x+1)**2 + sp.log(1-x)/(x-1) + sp.log(2)/x + 3
es = e.subs({sp.log(x): y1, sp.log(1-x): y2, sp.log(1+x): y3, sp.log(2): y4})
p = sp.Poly(es, y1, y2, y3, y4, domain=sp.QQ.frac_field(x))
print(p.terms())
for monom, coeff in p.terms():
    print(monom, "->", coeff, " zero?", sp.cancel(sp.sympify(coeff)) == 0)
