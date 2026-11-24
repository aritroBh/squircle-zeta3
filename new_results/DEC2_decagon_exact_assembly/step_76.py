import sympy as sp
x = sp.Symbol('x', positive=True)
Li2 = lambda z: sp.polylog(2,z)
L1 = Li2(x)+Li2(x/(x-1)) + sp.log(1-x)**2/2
d = sp.expand(sp.diff(L1,x))
d = d.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
# 1 - x/(x-1) = 1/(1-x) exactly -> log = -log(1-x)  (positive for 0<x<1)
d = d.xreplace({sp.log(1 - x/(x-1)): -sp.log(1-x)})
print("Landen d/dx final:", sp.simplify(d))
# also double-check the log argument identity exactly:
print("arg identity:", sp.simplify((1 - x/(x-1)) - 1/(1-x)))
