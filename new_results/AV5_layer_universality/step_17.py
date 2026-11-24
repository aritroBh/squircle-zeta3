import sympy as sp, mpmath, pickle
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, zeta, log, exp, quad, tanh, nstr
s = sp.symbols('s', positive=True)
f = {k: sp.sympify(v) for k,v in pickle.load(open('ckpt_fk_exprs.pkl','rb')).items()}
f3 = f[3]
print("f3 =", f3)
g1 = lambda x: log(1+exp(-2*x))
g3f = lambda x: mpf(2)/3*x**3*(tanh(x)-1)
pieces = lambda x: g3f(x) + g1(x)**3/6 + x*g1(x)**2/2 - x**2*g1(x)/2
func = sp.lambdify(s, f3, 'mpmath')
for x in ['0.3','1.0','2.0','5.0']:
    xx = mpf(x)
    print(x, nstr(func(xx),30), nstr(pieces(xx),30), nstr(func(xx)-pieces(xx),5))
