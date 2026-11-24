
# step_02c.py -- five-term relation, Suslin form
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')

def B(z):   # Bloch bracket {z} = delta11(S(Li_2(z))) = z ^ (1-z)
    return delta11(sym_G((0,1), z))

forms = {
 "Suslin: {x}-{y}+{y/x}-{(1-y)/(1-x)}+{x(1-y)/(y(1-x))}":
   [(x,1),(y,-1),(y/x,1),((1-y)/(1-x),-1),(x*(1-y)/(y*(1-x)),1)],
 "Zagier: {x}+{y}+{(1-x)/(1-xy)}+{(1-y)/(1-xy)}-{1-xy}":
   [(x,1),(y,1),((1-x)/(1-x*y),1),((1-y)/(1-x*y),1),(1-x*y,-1)],
 "v2: {x}+{y}-{(1-x)/(1-xy)}-{(1-y)/(1-xy)}+{1-xy}":
   [(x,1),(y,1),((1-x)/(1-x*y),-1),((1-y)/(1-x*y),-1),(1-x*y,1)],
}
for name, terms in forms.items():
    tot = {}
    for z, s in terms:
        tot = sym_add(tot, sym_scale(B(z), s))
    print(("OK  " if tot == {} else "no  "), name, "" if tot=={} else tot)
