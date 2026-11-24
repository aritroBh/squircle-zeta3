import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t10 = (-1-i)/2; t20 = (-1+i)/2
t11 = sp.Rational(3,4)*(1-i); t21 = sp.Rational(-3,4)-sp.Rational(3,4)*i
t12 = 1+sp.Rational(3,4)*i; t22 = 1-sp.Rational(3,4)*i
def rfun(dd, m):
    return 2*(dd-2*m) - 2*dd*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def rfunc(dd, m):  # conjugate (i -> -i)
    return 2*(dd-2*m) + 2*dd*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m):
    return sp.Rational(1,2)*( (3*dd-3)/2 - i*rfun(dd,m)/2 - (dd-2)**2 )
def R2c(dd, m):
    return sp.Rational(1,2)*( (3*dd-3)/2 + i*rfunc(dd,m)/2 - (dd-2)**2 )
d = 2*p-1
LHS = R2(d,p) - 1/(2*d-3)
rA0 = i*(2*d-1)/(p+1)
rB0 = (2*d-1)*(2*d+1)/((p+1)*p)
rB1 = -i*d*rB0
rB2 = R2(d+2,p+1)*rB0
# DL1 feed from (p+1,p-1) off-diag, depths 0..2
feed_DL1 = (p+1)*( t12*rA0/(2*d-1) + t11*(-i*(d-1))*rA0/(2*d-2) + t10*(R2(d+1,p+1)*rA0)/(2*d-3) )
# diagonal converted values
ed0 = t10*(4*p-3)/p
ed1 = (p+1)*( t10*rB1/(4*p-2) + t11*rB0/(4*p-1) )
ehat20 = (p+1)*( t10*rB2/(4*p-3) + t11*rB1/(4*p-2) + t12*rB0/(4*p-1) )
ehat21 = (p+1)*t10*(-rB0/(2*d+1))/(4*p-3)
ed2 = ehat20 - ehat21
feed_DL2 = p*( t22*ed0/(2*d-1) + t21*ed1/(2*d-2) + t20*ed2/(2*d-3) )
# fold at delta=2: m = 4p-6
ed0c = sp.conjugate(ed0)
ed1c = (p+1)*( t20*sp.conjugate(rB1)/(4*p-2) + (-t21)*sp.conjugate(rB0)/(4*p-1) )
ehat20c = (p+1)*( t20*(R2c(d+2,p+1)*rB0)/(4*p-3) + (-t21)*sp.conjugate(rB1)/(4*p-2) + t22*sp.conjugate(rB0)/(4*p-1) )
ehat21c = (p+1)*t20*(-rB0/(2*d+1))/(4*p-3)
ed2c = ehat20c - ehat21c
m = 4*p-6
g_m  = -( ed0c*t22/(4*p-3) + ed1c*(-t21)/(4*p-4) + ed2c*t20/(4*p-5) )
g_m1 = -( ed0c*(-t21)/(4*p-3) + ed1c*t20/(4*p-4) )
g_m2 = -( ed0c*t20/(4*p-3) )
feed_fold = p*( g_m - i*(4*p-7)*g_m1 - sp.binomial(4*p-6,2)*g_m2 )
RHS = feed_DL1 + feed_DL2 + feed_fold
print('simplifying E2b-adjacent...')
rem = sp.cancel(sp.together(LHS - RHS))
print('remainder:', rem)
