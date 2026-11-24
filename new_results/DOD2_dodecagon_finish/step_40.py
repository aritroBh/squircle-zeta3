import sympy as sp, inspect
print("sympy version:", sp.__version__)
from sympy.matrices import reductions
src = inspect.getsource(sp.Matrix.lll)
print(src[:3000])
