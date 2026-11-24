import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t10 = (-1-i)/2; t20 = (-1+i)/2; t11 = sp.Rational(3,4)*(1-i); t21 = sp.Rational(-3,4)-sp.Rational(3,4)*i
d = 2*p-1
LHS = -i*(d-2)
rA0 = i*(2*d-1)/(p+1)
rA1 = -i*(d-1)*rA0
rB0 = (2*d-1)*(2*d+1)/((p+1)*p)     # rho_{p+1,p}(0)/rho_{p,p-1}(0)  (real)
rB1 = -i*d*rB0
ehat00 = t10*(4*p-3)/p
ehat10 = (p+1)*( t10*rB1/(4*p-2) + t11*rB0/(4*p-1) )
feed_DL1 = (p+1)*( t10*rA1/(2*d-2) + t11*rA0/(2*d-1) )
feed_DL2 = p*( t20*ehat10/(2*d-2) + t21*ehat00/(2*d-1) )
ehat00c = sp.conjugate(ehat00)
ehat10c = (p+1)*( t20*sp.conjugate(rB1)/(4*p-2) + t21*sp.conjugate(rB0)/(4*p-1) )
g_m1c = -( ehat00c*t21/(4*p-1) + ehat10c*t20/(4*p-2) )
g_m0c = -( ehat00c*t20/(4*p-1) )
feed_fold = p*( g_m1c + (4*p-6)*g_m0c )
RHS = feed_DL1 + feed_DL2 + feed_fold
print('E1-adjacent simplified:', sp.cancel(sp.together(LHS - RHS)))
