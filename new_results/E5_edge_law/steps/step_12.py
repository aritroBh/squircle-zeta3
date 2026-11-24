# step_12.py: per-edge F_k(L,h,s0) = int_0^1 E_k dt, k=2..5, for a PARALLELOGRAM edge
# (alpha=beta=1, mu0=mu1=0, l0=X, l1=Y), exact in {1,z2,z3,z2^2,z2z3,z5}, symbolic in L,h,s0.
# Then P^(k)(1) = k! * sum over 4 edges F_k.  Verify: F_2 = z2 h^2/L, F_3 = h^2/L (z3 - 2 z2)  [E1].
import sympy as sp, pickle, sys
sys.path.insert(0, 'steps')
from step_03 import J, z2, z3, z5

t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1, X, Y = sp.symbols('l0 l1 X Y')
R = pickle.load(open('steps/step_01_reduced.pkl', 'rb'))
loc = {'t': t, 'alpha': al, 'beta': be, 'L': L, 'h': h, 's0': s0, 'l0': l0, 'l1': l1}
E = [sp.sympify(s, locals=loc) for s in R['E']]
F = []
for k in range(6):
    ex = sp.expand(E[k].subs({al: 1, be: 1, l0: X, l1: Y}))
    P = sp.Poly(ex, t, X, Y)
    tot = 0
    for (n, a, b), c in P.terms():
        tot += c*J(n, 0, a, b)
    F.append(sp.expand(tot))
    print(f"F_{k} =", sp.collect(F[k], [z5, z2*z3, z3, z2]))
# checks vs E1
print("F_2 - z2 h^2/L =", sp.simplify(F[2] - z2*h**2/L))
print("F_3 - h^2/L (z3-2z2) =", sp.simplify(F[3] - h**2/L*(z3 - 2*z2)))
# structure of F_4, F_5: split by h-power
for k in [4, 5]:
    expr = F[k]
    deg = sp.Poly(sp.expand(expr*L**6), L)  # just to see L powers
    print(f"F_{k} as poly in L:", sp.Poly(sp.expand(expr), L).terms() if False else "")
    c = sp.collect(sp.expand(expr), [h], evaluate=False)
    print(f"F_{k} by h-power:", {str(key): sp.collect(val, [z5, z2*z3, z3, z2]) for key, val in c.items()})
pickle.dump(dict(F=[str(x) for x in F]), open('steps/step_12_F.pkl', 'wb'))
with open('steps/step_12_out.txt', 'w') as f:
    for k in range(6):
        f.write(f"F_{k} = {sp.collect(F[k], [z5, z2*z3, z3, z2])}\n")
print("saved steps/step_12_F.pkl / step_12_out.txt")
