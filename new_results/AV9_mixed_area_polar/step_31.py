
import sympy as sp
z2 = sp.zeta(2)
d = sp.Symbol('piD')
e = sp.pi**4/720 + sp.pi**2/2
e2 = sp.expand(e.subs(sp.pi, d))
print("dummy subs:", e2.subs(d**2, 6*z2))
# So conversion: e -> subs(pi,d) -> subs(d**2,6*z2)
