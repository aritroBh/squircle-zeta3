import sympy as sp, time
u = sp.symbols('u', positive=True)
s = sp.sin(u); c = sp.cos(u); s2 = s**2; c2 = 1-s2
ls = sp.log(s); lc2 = sp.log(c2)
A = 2*s2*ls - (s2-1)*lc2
H1u = A/(4*c2**sp.Rational(3,2))*c
t0=time.time()
# quick sanity on H1 (known closed form)
i1 = sp.integrate(H1u, u)
print("H1 int:", time.time()-t0, "s")
print(sp.simplify(i1))
