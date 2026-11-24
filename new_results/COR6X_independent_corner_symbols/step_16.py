dbg = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')
z1 = (1-x)/(1-x*y); z2 = (1-y)/(1-x*y); z3 = 1-x*y
for name, z in [("x",x),("y",y),("z1",z1),("z2",z2),("1-xy",z3)]:
    SG = sym_G((0,1), z)
    print(name, ": S(G(0,1;z)) =", SG)
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",dbg],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
