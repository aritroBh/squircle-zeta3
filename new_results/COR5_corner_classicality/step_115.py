import sympy as sp
print(sp.__version__)
import sympy.polys.numberfields as nf
print([x for x in dir(nf) if 'inim' in x.lower()])
