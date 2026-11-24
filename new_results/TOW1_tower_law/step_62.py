import sympy as sp
d, p = sp.symbols('d p', integer=True, positive=True)
i = sp.I
t10 = (-1-i)/2; t20 = (-1+i)/2
t11 = sp.Rational(3,4)*(1-i); t21 = sp.Rational(-3,4)-sp.Rational(3,4)*i
t12 = 1+sp.Rational(3,4)*i; t22 = 1-sp.Rational(3,4)*i
def rfun(dd, m):
    return 2*(dd-2*m) - 2*dd*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m):
    return sp.Rational(1,2)*( (3*dd-3)/2 - i*rfun(dd,m)/2 - (dd-2)**2 )
LHS = R2(d,p) - 1/(2*d-3)
rA0 = i*(2*d-1)/(p+1)
rB0 = -i*(2*d-1)/(d-p+1)
T0  = (p+1)*t12*rA0/(2*d-1)
T1  = (p+1)*t11*(-i*(d-1))*rA0/(2*d-2)
T2  = (p+1)*t10*(R2(d+1,p+1)*rA0)/(2*d-3)
T0b = (d-p+1)*t22*rB0/(2*d-1)
T1b = (d-p+1)*t21*(-i*(d-1))*rB0/(2*d-2)
T2b = (d-p+1)*t20*(R2(d+1,p)*rB0)/(2*d-3)
RHS = T0+T1+T2+T0b+T1b+T2b
print('simplifying E2b non-adjacent...')
rem = sp.cancel(sp.together(LHS - RHS))
print('remainder:', rem)
