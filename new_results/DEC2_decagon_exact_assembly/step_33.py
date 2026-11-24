import sympy as sp
# exact trig values at i*pi/5
for i in range(1,6):
    c = sp.nsimplify(sp.cos(i*sp.pi/5))
    s = sp.nsimplify(sp.sin(i*sp.pi/5))
    print(i, "cos =", sp.simplify(sp.cos(i*sp.pi/5)), " sin =", sp.simplify(sp.sin(i*sp.pi/5)))
