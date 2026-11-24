import sympy as sp
p = sp.symbols('p', positive=True)
pip = 2*sp.pi/(p*sp.sin(sp.pi/p))
N = 4/(p*pip)*sp.beta(3/p, 1-1/p)
N2 = sp.diff(N, p, 2).subs(p, 2)
# manually expand the big log
N2e = sp.expand(N2)
N2e = N2e.replace(lambda e: isinstance(e, sp.log) and e.args[0].is_Integer,
                  lambda e: sum(m*sp.log(q) for q,m in sp.factorint(e.args[0]).items()))
N2e = sp.simplify(N2e)
print("N''(2) =", N2e)
print(sp.N(N2e,30))
