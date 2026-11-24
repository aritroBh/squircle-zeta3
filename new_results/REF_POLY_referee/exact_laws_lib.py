import sympy as sp
lt, l1 = sp.symbols('lt l1')
xx, yy = sp.symbols('xx yy', positive=True)
Bfun = sp.gamma(xx)*sp.gamma(yy)/sp.gamma(xx+yy)
_tab = {}
def mom(n, a, b):
    key = (n, a, b)
    if key not in _tab:
        v = sp.diff(Bfun, xx, a, yy, b).subs({xx: n+1, yy: 1})
        _tab[key] = sp.expand(sp.simplify(sp.expand_func(v)))
    return _tab[key]
