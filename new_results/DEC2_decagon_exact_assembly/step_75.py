import sympy as sp
x = sp.Symbol('x', positive=True)
Li2 = lambda z: sp.polylog(2,z)
L1 = Li2(x)+Li2(x/(x-1)) + sp.log(1-x)**2/2
d = sp.expand(sp.diff(L1,x))
d = d.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
print("atoms:", d.atoms(sp.log))
key = sp.log(-1/(x-1))
d2 = d.xreplace({key: -sp.log(1-x)})
print("after xreplace:", sp.simplify(d2))
# alternative: also try the form sympy might use
print("check key in d:", d.has(key))
# maybe it's log(1/(1-x))
for cand in [sp.log(1/(1-x)), sp.log(1/(1-x))]:
    print(cand, d.has(cand))
