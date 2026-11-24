import sympy as sp
t, z = sp.symbols('t z')
f = {1: sp.Integer(1), 2: t, 3: 1+3*t, 4: 2*t+3*t**2, 5: 2+15*t+15*t**2, 6: 17*t+60*t**2+45*t**3, 7: 17+231*t+525*t**2+315*t**3}
a1,a0,b1,b0,c1,c0 = sp.symbols('a1 a0 b1 b0 c1 c0')
eqs = []
for d in [4,5,6,7]:
    expr = sp.expand(f[d] - ((a1*t+a0)*f[d-1] + (b1*t+b0)*f[d-2] + (c1*t+c0)*f[d-3]))
    for cf in sp.Poly(expr, t).all_coeffs():
        eqs.append(sp.Eq(cf, 0))
sol = sp.solve(eqs, [a1,a0,b1,b0,c1,c0], dict=True)
print("order-3 sol:", sol)
a2,a1_,a0_,b2,b1_,b0_ = sp.symbols('a2 a1 a0 b2 b1 b0')
eqs2 = []
for d in [4,5,6,7]:
    expr = sp.expand(f[d] - ((a2*t**2+a1_*t+a0_)*f[d-1] + (b2*t**2+b1_*t+b0_)*f[d-2]))
    for cf in sp.Poly(expr, t).all_coeffs():
        eqs2.append(sp.Eq(cf, 0))
sol2 = sp.solve(eqs2, [a2,a1_,a0_,b2,b1_,b0_], dict=True)
print("order-2 quad sol:", sol2)
