import sympy as sp
from sympy.functions.combinatorial.numbers import stirling
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
bad = []
for r in range(1,5):
    for n in range(0, 25):
        sv = stirling(n, n-r, kind=2)
        pv = P[r].subs(x, n)
        if sp.Integer(sv) != pv:
            bad.append((r, n, sv, pv))
print(len(bad))
for b in bad[:20]: print(b)
