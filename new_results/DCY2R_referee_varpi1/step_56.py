import sympy as sp
t = sp.symbols('t')
h1,h2,h3,l3,l4,l5 = sp.symbols('h1 h2 h3 l3 l4 l5')
PH = 1 + h1*t + h2*t**2 + h3*t**3
S = l3*t**3 + l4*t**4 + l5*t**5
QS = 1 + S + S**2/2 + S**3/6
full = sp.expand(PH*QS)
# their kept P
P = (1 + h1*t + h2*t**2 + (h3+l3)*t**3 + (l4+h1*l3)*t**4 + (l5+h1*l4+h2*l3)*t**5
     + (l3**2/2)*t**6 + (l3*l4+h1*l3**2/2)*t**7 + (l3**3/6)*t**9)
# their drop list
drop = {6: h1*l5+h2*l4+h3*l3, 7: h2*l5+h3*l4,
        8: h1*l3*l4+h2*l3**2/2+h3*l5, 9: h2*l3*l4+h3*l3**2/2,
        10: h1*l3**3/6, 11: h2*l3**3/6, 12: h3*l3**3/6}
rest = sp.expand(full - P)
print("monomials of full-P:", sorted(rest.as_poly(t).monoms()))
ok = True
for p,c in drop.items():
    coeff = rest.coeff(t, p)
    if sp.simplify(coeff - c) != 0:
        ok = False; print(f"t^{p}: drop coeff mismatch: {coeff} vs {c}")
# remaining monomials not in drop and not covered by beta2/beta3 claims:
# from QS - kept: (l3l5+l4^2/2) t^8 + l4l5 t^9 + l5^2 t^10/2  [S^2/2 leftovers]
# and S^3/6 leftovers: 3 l3^2 l4 t^10/6 ... etc; and PH * those
resid_terms = {}
for (p,), c in rest.as_poly(t).terms():
    resid_terms[p] = c
print("residual (full - P - droplist):")
for p in sorted(resid_terms):
    c = sp.simplify(resid_terms[p] - drop.get(p,0))
    if c != 0:
        print(f"  t^{p}: {c}")
