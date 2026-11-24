import sympy as sp
x = sp.Symbol('x', positive=True)
Li2 = lambda z: sp.polylog(2,z)
L1 = Li2(x)+Li2(x/(x-1)) + sp.log(1-x)**2/2
d = sp.expand(sp.diff(L1,x))
d = d.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
logs = set()
d.find(lambda f: isinstance(f, sp.log), lambda f: logs.add(f))
print("logs present:", logs)
# 1 - x/(x-1) = 1/(1-x); replace accordingly
rep = {}
for f in logs:
    arg = f.args[0]
    arg2 = sp.simplify(arg)
    if sp.simplify(arg2 - 1/(1-x)) == 0:
        rep[f] = -sp.log(1-x)
    elif sp.simplify(arg2 - (1-x)) == 0:
        rep[f] = sp.log(1-x)
print("replacement:", rep)
d2 = sp.simplify(d.subs(rep, simultaneous=True))
print("Landen d/dx final:", d2)
