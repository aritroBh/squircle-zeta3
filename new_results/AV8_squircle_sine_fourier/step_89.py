import sympy as sp
p = sp.symbols('p', positive=True)
pip = 2*sp.pi/(p*sp.sin(sp.pi/p))
N = 4/(p*pip)*sp.beta(3/p, 1-1/p)
N3raw = sp.diff(N, p, 3).subs(p, 2)
N3e = N3raw.replace(lambda e: isinstance(e, sp.log) and e.args[0].is_Integer,
                lambda e: sum(mm*sp.log(q) for q,mm in sp.factorint(e.args[0]).items()))
N3e = N3e.subs({sp.polygamma(1, sp.Rational(3,2)): sp.pi**2/2-4,
                sp.polygamma(2, sp.Rational(3,2)): -14*sp.zeta(3)+16})
N3e = sp.simplify(N3e)
print("N'''(2) =", N3e)
print(sp.N(N3e,30))
