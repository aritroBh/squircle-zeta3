import sympy as sp
p = sp.symbols('p', positive=True)
pip = 2*sp.pi/(p*sp.sin(sp.pi/p))
N = 4/(p*pip)*sp.beta(3/p, 1-1/p)
N3raw = sp.diff(N, p, 3).subs(p, 2)
print("numeric raw:", sp.N(N3raw, 30))
# expand big logs and substitute polygamma values exactly
N3e = N3raw.replace(lambda e: isinstance(e, sp.log) and e.args[0].is_Integer,
                lambda e: sum(mm*sp.log(q) for q,mm in sp.factorint(e.args[0]).items()))
print("polygamma args present:", set(sp.preorder_traversal(N3e)) and [a for a in sp.preorder_traversal(N3e) if isinstance(a, sp.polygamma)])
