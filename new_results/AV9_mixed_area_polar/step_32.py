
import sympy as sp
z2 = sp.zeta(2); d = sp.Symbol('piD')
e = sp.pi**4/720 + sp.pi**2/2
e2 = e.subs(sp.pi, d)
print("e2 =", e2, type(e2))
print("free:", e2.free_symbols)
e3 = e2.subs(d**2, 6*z2)
print("e3 =", e3)
print("e3 direct:", (d**4/720).subs(d**2, 6*z2))
