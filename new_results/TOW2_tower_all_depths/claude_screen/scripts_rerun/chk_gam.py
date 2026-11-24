import sympy as sp
d = sp.Symbol('d'); s = sp.Symbol('s', positive=True, integer=True); h = sp.Rational(1,2)
G = sp.gamma
ratio = (G(d+h)/G(d-h)) * (G(d-s-h)/G(d-s+h))
print('Gamma ratio simplifies (via recurrence Gamma(x+1)=x Gamma(x)):', sp.simplify(sp.gammasimp(sp.expand_func(G(d+h)) / G(d-h)) - (d-h)), sp.simplify(sp.expand_func(G(d-s+h))/G(d-s-h)))
print('=> ratio =', sp.simplify((d-h)/(d-s-h)), ' equals (2d-1)/(2d-2s-1):', sp.simplify((d-h)/(d-s-h) - (2*d-1)/(2*d-2*s-1)) == 0)
