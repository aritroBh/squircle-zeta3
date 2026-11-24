import sympy as sp, tow2lib as T
d, m = T.d, T.m; i = sp.I
r = 2*(d-2*m) - 2*d*i - 4*(d-1)*(d-2*m)/(2*d-3)
a2 = sp.Rational(1,2)*((3*d-3)/sp.Integer(2) - i*r/2 - (d-2)**2)
print('arho_simple(2) == TOW1 rho(2)/rho(0):', sp.cancel(sp.together(T.arho_simple(2) - a2)) == 0)
print('arho_simple(1) == -i(d-2):', sp.cancel(T.arho_simple(1) + i*(d-2)) == 0)
