
# step_02b.py -- five-term relation: find the exact relation among the five delta11 vectors
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')

args = [x, y, (1-x)/(1-x*y), (1-y)/(1-x*y), 1-x*y]
vecs = []
keys = []
for z in args:
    v = delta11(sym_G((0,1), z))   # S(Li_2(z)) = -[1-z][z] -> delta11 = z^(1-z)
    vecs.append(v)
    for k in v:
        if k not in keys: keys.append(k)
print("wedge basis size:", len(keys))
import sympy as sp
M = sp.Matrix([[v.get(k, 0) for k in keys] for v in vecs])
print("rank:", M.rank())
ns = M.nullspace()
print("nullspace:", ns)
for rel in ns:
    combo = [str(args[i]) for i in range(5) if rel[i] != 0]
    print("relation coefficients:", list(rel), "on args:", combo)

# classical five-term: [x] + [y] + [(1-x)/(1-xy)] + [(1-y)/(1-xy)] + [1-xy] = 0 ?
# with [z] meaning delta11(S(Li_2(z))): test all-plus
tot = {}
for v in vecs: tot = sym_add(tot, v)
print("all-plus sum == 0:", tot == {})
