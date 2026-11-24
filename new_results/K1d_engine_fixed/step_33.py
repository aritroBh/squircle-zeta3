import sympy as sp
t, k1, k2 = sp.symbols('t k1 k2', positive=True)
# Prototype: the new (c1,c2)=(1,1) weight-2 integral
I = sp.integrate(sp.log(1+k1*t)*sp.log(1+k2*t), (t, 0, 1))
print("int_0^1 log(1+k1 t) log(1+k2 t) dt =")
print(sp.simplify(I))
