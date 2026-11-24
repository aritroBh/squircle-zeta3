import sympy as sp
p, y = sp.symbols('p y', positive=True)
# g(p) = -(1/p) ln(1-y^p); H = exp(g); H_k = d^k H/dp^k at p=2
g = -sp.log(1-y**p)/p
g1 = sp.diff(g,p); g2 = sp.diff(g1,p); g3 = sp.diff(g2,p)
H = sp.exp(g)
H1 = sp.diff(H,p); H2 = sp.diff(H1,p); H3 = sp.diff(H2,p)
def at2(expr):
    return sp.simplify(expr.subs(p,2))
g1_2 = at2(g1); g2_2 = at2(g2); g3_2 = at2(g3)
print("g1|2 =", g1_2)
print("g2|2 =", sp.simplify(g2_2))
print("g3|2 =", sp.simplify(g3_2))
H1_2 = sp.simplify(at2(H1)); H2_2 = sp.simplify(at2(H2)); H3_2 = sp.simplify(at2(H3))
print("H1|2 =", H1_2)
print("H2|2 =", H2_2)
print("H3|2 =", H3_2)
