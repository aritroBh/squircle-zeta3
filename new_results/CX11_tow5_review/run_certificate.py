"""Run local certificates with exact Gaussian-rational cancellation.
No identity is changed. Explicit coefficient domains avoid SymPy's large EX expansions.
"""
import runpy, sys
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.domains import QQ_I
orig_cancel=sp.cancel
# All cancellation inputs in these certificates are rational over Q(i).
def gaussian_cancel(expr,*args,**kwargs):
    if args or kwargs: return orig_cancel(expr,*args,**kwargs)
    syms=sorted(expr.free_symbols,key=str)
    if not syms: return sp.simplify(expr)
    K,*_=field(syms,QQ_I)
    return K.from_expr(expr).as_expr()
sp.cancel=gaussian_cancel
script=sys.argv[1];sys.argv=sys.argv[1:]
if script=='verify_L.py':
    import tow5lib as T5, tow2lib as T2
    for l in T2.sl: T2.sl[l]=sp.cancel(T2.sl[l])
runpy.run_path(script,run_name='__main__')
