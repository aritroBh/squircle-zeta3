import sympy as sp
w, ca = sp.symbols('w ca')
Q = sp.sqrt(1 - 2*ca*w + w**2)
# (i) symmetry: q(1/w) = w^2 q(w) for q = (1+w)/Q^3
q = (1+w)/Q**3
print("q(1/w) - w^2 q(w) =", sp.simplify(q.subs(w, 1/w) - w**2*q))
# h(1/w) = w^2 h(w) for h = [(1+w)/Q - 1]/w
h = ((1+w)/Q - 1)/w
print("h(1/w) - w^2 h(w) =", sp.simplify(h.subs(w, 1/w) - w**2*h))
# (ii) primitive: d/dw[1-(1+w)/Q] = (1+ca)(w-1)/Q^3
expr = sp.simplify(sp.diff(1-(1+w)/Q, w) - (1+ca)*(w-1)/Q**3)
print("diff[1-(1+w)/Q] - (1+ca)(w-1)/Q^3 =", sp.radsimp(expr))
# (v) V(w)=(w-1)/((1-ca)Q) primitive of (1+w)/Q^3
print("diff V - (1+w)/Q^3 =", sp.radsimp(sp.diff((w-1)/((1-ca)*Q), w) - (1+w)/Q**3))
# Q(1/w) = Q(w)/w
print("Q(1/w) - Q(w)/w =", sp.simplify(Q.subs(w,1/w) - Q/w))
