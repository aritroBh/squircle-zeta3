import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, gamma, pi
def F(y):
    return quad(lambda w: w**(-y-1)*(1+w)*(1/sqrt(1+w**2)-1), [0,1])
def E(y):
    return gamma(1+y/2)*gamma(mpf('0.5')-y/2)/sqrt(pi)
def Beta_piece(y):  # (E(y)-E(-y))/(2y) via Beta integrals directly
    return (E(y)-E(-y))/(2*y)
def beta_split(y):
    I1 = quad(lambda t: t**(y/2)*(1-t)**(-mpf('0.5')-y/2), [0, mpf('0.5'), 1])
    I2 = quad(lambda t: t**(-y/2)*(1-t)**(-mpf('0.5')+y/2), [0, mpf('0.5'), 1])
    return (I1-I2)/(4*y)
def w_form(y):  # claimed: (1/(2y))∫(w^y-w^{-y})(w-1)/(1+w^2)^{3/2} dw
    return quad(lambda w: (w**y - w**(-y))*(w-1)/(1+w**2)**mpf('1.5'), [0,1])/(2*y)
for y in [mpf('0.3'), mpf('0.55')]:
    print('y=',y)
    print('  F+F-/2      =', mpmath.nstr((F(y)+F(-y))/2, 30))
    print('  (E-E-)/(2y) =', mpmath.nstr(Beta_piece(y), 30))
    print('  beta_split  =', mpmath.nstr(beta_split(y), 30))
    print('  w_form      =', mpmath.nstr(w_form(y), 30))
