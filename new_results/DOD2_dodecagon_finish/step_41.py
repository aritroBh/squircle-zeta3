import inspect
from sympy.polys.matrices import domainmatrix
from sympy.polys.matrices.lll import ddm_lll, ddm_lll_transform
src = inspect.getsource(ddm_lll)
print(src[:4000])
