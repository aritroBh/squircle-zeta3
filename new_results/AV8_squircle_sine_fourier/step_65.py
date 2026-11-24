import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
H = sp.harmonic
d = (-m**2*H(m) + m**2*H(2*m) + m**2*H(m-1) - m**2*H(2*m-2) + m*H(m)/2 - m*H(2*m)/2 - m*H(m-1)/2 + m*H(2*m-2)/2 + sp.Rational(1,4))/(m*(2*m-1))
d2 = d.subs({H(m-1): H(m)-1/m, H(2*m-2): H(2*m)-1/(2*m-1)-1/(2*m)})
print("difference after shifts:", sp.simplify(d2))
