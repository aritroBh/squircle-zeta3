import sympy as sp
p, y = sp.symbols('p y', positive=True)
g = -sp.log(1-y**p)/p
H = sp.exp(g)
Hs = [sp.simplify(sp.diff(H,p,k).subs(p,2)) for k in [1,2,3]]
s2, c2, ls, lc2 = sp.symbols('s2 c2 ls lc2')
rep = {y**2: s2, (1-y**2): c2, sp.log(y): ls, sp.log(1-y**2): lc2}
lines = []
for k,Hk in enumerate(Hs,1):
    num = sp.fraction(Hk)[0]*1
    den = sp.fraction(Hk)[1]
    expr = sp.simplify(Hk*(1-y**2)**(sp.Rational(2*k+1,2)))
    expr = sp.expand(expr)
    expr2 = expr.subs({y**6: s2**3, y**4: s2**2, y**2: s2, sp.log(y): ls, sp.log(1-y**2): lc2})
    expr2 = expr2.subs({(1-y**2): c2})
    expr2 = sp.expand(expr2).subs({y**2: s2, 1-y**2: c2, sp.log(y): ls, sp.log(1-y**2): lc2})
    print(k, ":", expr2)
    lines.append(str(expr2))
# check power of denominator
for k,Hk in enumerate(Hs,1):
    print(k, sp.denom(Hk))
