import sympy as sp
x = sp.symbols('x')
# How does sympy differentiate polylog?
e = sp.polylog(3, (1-x)/(1+x))
print("diff Li3:", sp.diff(e, x))
e2 = sp.polylog(2, 2*x/(1+x))
print("diff Li2:", sp.diff(e2, x))
print("polylog(1,z):", sp.polylog(1, x), "-> auto?", sp.simplify(sp.polylog(1, x)))
