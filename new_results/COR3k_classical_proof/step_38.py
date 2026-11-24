import sympy as sp
x = sp.symbols('x', positive=True)
p1 = sp.polylog(2, (1-x)/(1+x))
print("constructed arg:", p1.args[1], "| hash:", hash(p1.args[1]))
# what dpoly produces:
d = sp.diff(sp.polylog(3, (1-x)/(1+x)), x)
pl = [t for t in sp.preorder_traversal(d) if isinstance(t, sp.polylog)][0]
print("from diff arg:  ", pl.args[1], "| hash:", hash(pl.args[1]))
print("structurally equal:", p1.args[1] == pl.args[1], "| p1 == pl:", p1 == pl)
e = sp.expand((1-x)/(1+x))
print("expand((1-x)/(1+x)) =", e, "| cancel:", sp.cancel((1-x)/(1+x)))
print("together:", sp.together(pl.args[1]))
