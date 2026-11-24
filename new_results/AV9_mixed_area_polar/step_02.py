
# step_02.py : derive h+h'' for h=(c^q+s^q)^{1/q}; build bivariate W(q,u) formula; validate
import sympy as sp
th, q = sp.symbols('th q', positive=True)
c, s = sp.cos(th), sp.sin(th)
A = c**q + s**q
h = A**(1/q)
R = sp.simplify(h + sp.diff(h, th, 2))
# manual: compute log derivative
# h'/h = (1/q) A'/A ; A' = q(c s^{q-1} - s c^{q-1})
Ap = q*(c*s**(q-1) - s*c**(q-1))
App = sp.diff(Ap, th)
# h''/h = (1/q)(App/A - (Ap/A)^2) + ((1/q)(Ap/A))^2 = (1/q)App/A + (1/q)(1/q-1)(Ap/A)^2
Rmanual = h*(1 + sp.Rational(1,1)/q*App/A + (1/q)*(1/q-1)*(Ap/A)**2)
# simplify App
App_s = sp.simplify(App)
print("A'' =", sp.trigsimp(App_s))
# guess closed form: R = (q-1) c^{q-2} s^{q-2} A^{1/q-2} * (something)
guess = (q-1)*c**(q-2)*s**(q-2)*A**(1/q-2)
ratio = sp.simplify(Rmanual/guess)
print("Rmanual / guess =", sp.trigsimp(sp.expand(ratio)))
