import inspect
from sympy.polys.matrices import lll as lllmod
src = inspect.getsource(lllmod._ddm_lll)
print(src)
