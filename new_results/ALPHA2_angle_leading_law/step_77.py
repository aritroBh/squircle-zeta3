import sympy as sp
x = sp.symbols('x')
P = {1: sp.Integer(1), 2: 1-x, 3: 2-x, 4: 5-6*x+x**2, 5: 16-13*x+x**2,
     6: 61-89*x+29*x**2-x**3, 7: 272-297*x+60*x**2-x**3, 8: 1385-2396*x+1134*x**2-124*x**3+x**4}
a,b,c,e,g,h = sp.symbols('a b c e g h')
eqs = []
for d in [5,6,7,8]:
    expr = sp.expand(P[d] - ((a+b*x)*P[d-1] + (c+e*x)*P[d-2] + (g+h*x)*P[d-3]))
    for cf in sp.Poly(expr, x).all_coeffs():
        eqs.append(sp.Eq(cf, 0))
sol = sp.solve(eqs, [a,b,c,e,g,h], dict=True)
print("P order-3:", sol)
# also try order-2 with quadratic coeff
a2,a1,a0,b2,b1,b0 = sp.symbols('a2 a1 a0 b2 b1 b0')
eqs2 = []
for d in [4,5,6,7,8]:
    expr = sp.expand(P[d] - ((a2*x**2+a1*x+a0)*P[d-1] + (b2*x**2+b1*x+b0)*P[d-2]))
    for cf in sp.Poly(expr, x).all_coeffs():
        eqs2.append(sp.Eq(cf, 0))
print("P order-2 quad:", sp.solve(eqs2, [a2,a1,a0,b2,b1,b0], dict=True))
# and order-4 linear
cs = sp.symbols('a1 a0 b1 b0 c1 c0 e1 e0')
a1,a0,b1,b0,c1,c0,e1,e0 = cs
eqs3 = []
for d in [6,7,8]:
    expr = sp.expand(P[d] - ((a1*x+a0)*P[d-1] + (b1*x+b0)*P[d-2] + (c1*x+c0)*P[d-3] + (e1*x+e0)*P[d-4]))
    for cf in sp.Poly(expr, x).all_coeffs():
        eqs3.append(sp.Eq(cf, 0))
print("P order-4 lin:", sp.solve(eqs3, cs, dict=True))
