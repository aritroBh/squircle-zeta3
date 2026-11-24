import sympy as sp
import mpmath
from mpmath import mp, mpf, sqrt, pi, zeta
mp.dps = 70

tau, s = sp.symbols('tau s', positive=True)
th = sp.tanh(s)
g = {
 1: sp.log(2*sp.cosh(s)) - s,
 3: sp.Rational(2,3)*s**3*(th-1),
 5: sp.Rational(2,3)*s**5*(th-1) + sp.Rational(2,9)*s**6*(1-th**2),
 7: (-sp.Rational(8,81)*s**9*th*(1-th**2) + sp.Rational(4,9)*s**8*(1-th**2)
     + sp.Rational(244,315)*s**7*(th-1)),
}
G = sum(tau**k*gk for k,gk in g.items())
eG = sp.series(sp.exp(G), tau, 0, 10).removeO().expand()
css = sp.series(sp.cos(tau*s)+sp.sin(tau*s), tau, 0, 10).removeO().expand()
prod = sp.expand(eG.subs(sp.exp(G), 1)*0 + sp.expand((eG-1)*css))
f = {}
for k in range(1,10):
    fk = sp.expand(prod.coeff(tau,k))
    f[k] = fk
    print(f"f_{k}(s) built, {len(sp.Add.make_args(fk)) if fk.is_Add else 1} terms")
sp.save = None
import pickle
pickle.dump({k: str(v) for k,v in f.items()}, open('ckpt_fk_exprs.pkl','wb'))
