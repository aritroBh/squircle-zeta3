import sympy as sp
t = sp.symbols('t')
f = {1: sp.Integer(1), 2: t, 3: 1+3*t, 4: 2*t+3*t**2, 5: 2+15*t+15*t**2,
     6: 17*t+60*t**2+45*t**3, 7: 17+231*t+525*t**2+315*t**3, 8: 62*t+378*t**2+630*t**3+315*t**4}
# order-4 recurrence with t-linear coefficients
cs = sp.symbols('a1 a0 b1 b0 c1 c0 e1 e0')
a1,a0,b1,b0,c1,c0,e1,e0 = cs
eqs = []
for d in [5,6,7,8]:
    expr = sp.expand(f[d] - ((a1*t+a0)*f[d-1] + (b1*t+b0)*f[d-2] + (c1*t+c0)*f[d-3] + (e1*t+e0)*f[d-4]))
    for cf in sp.Poly(expr, t).all_coeffs():
        eqs.append(sp.Eq(cf, 0))
sol = sp.solve(eqs, cs, dict=True)
print("order-4 t-linear:", sol)
# order-2 with cubic-in-t coefficients
cs2 = sp.symbols('a3 a2 a1 a0 b3 b2 b1 b0')
eqs2 = []
for d in [4,5,6,7,8]:
    a3,a2,a1,a0,b3,b2,b1,b0 = cs2
    expr = sp.expand(f[d] - ((a3*t**3+a2*t**2+a1*t+a0)*f[d-1] + (b3*t**3+b2*t**2+b1*t+b0)*f[d-2]))
    for cf in sp.Poly(expr, t).all_coeffs():
        eqs2.append(sp.Eq(cf, 0))
sol2 = sp.solve(eqs2, cs2, dict=True)
print("order-2 t-cubic:", sol2)
