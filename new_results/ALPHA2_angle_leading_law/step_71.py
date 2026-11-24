import sympy as sp
w = sp.symbols('w', positive=True); ca = sp.symbols('ca')
Q = sp.sqrt(1 - 2*ca*w + w**2)
q = (1+w)/Q**3
e1 = sp.radsimp(sp.simplify(q.subs(w, 1/w).rewrite(sp.Pow) - w**2*q))
print("q symmetry residual:", sp.simplify(e1))
h = ((1+w)/Q - 1)/w
e2 = sp.radsimp(sp.simplify(h.subs(w, 1/w) - w**2*h))
print("h symmetry residual:", sp.simplify(e2))
V = (w-1)/((1-ca)*Q)
e3 = sp.radsimp(sp.simplify(sp.diff(V, w) - (1+w)/Q**3))
print("V primitive residual:", sp.simplify(e3))
e4 = sp.simplify(Q.subs(w,1/w) - Q/w)
print("Q(1/w)-Q/w:", sp.radsimp(e4))
