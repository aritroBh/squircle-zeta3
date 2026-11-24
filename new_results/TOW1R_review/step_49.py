import sympy as sp
i=sp.I; sq2=sp.sqrt(2); a=1+sq2; b=sq2-1; v=sp.symbols('v')
f = 2*v/(v**2-a**2)+2*v/(v**2-b**2)
dfdv = sp.diff(f,v)
print("d/dv[2v/(v^2-a^2)+2v/(v^2-b^2)] at i =", sp.simplify(dfdv.subs(v,i)))
print("negated                                                    =", sp.simplify(-dfdv.subs(v,i)))
# So in (logA)'' the term is -d/dv[...]; what value?
# The proof text claims -1. Actual:
# also recompute via my closed form -2(v^2+c^2)/(v^2-c^2)^2
for c in (a,b):
    print(f"c={c}: d/dv[2v/(v^2-c^2)] at i direct =", sp.simplify(sp.diff(2*v/(v**2-c**2),v).subs(v,i)),
          "  my -2(v^2+c^2)/(v^2-c^2)^2 =", sp.simplify(-2*(i**2+c**2)/(i**2-c**2)**2))
