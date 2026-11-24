import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t10 = (-1-i)/2; t20 = (-1+i)/2; t11 = sp.Rational(3,4)-sp.Rational(3,4)*i; t21 = sp.Rational(-3,4)-sp.Rational(3,4)*i
# E1 adjacent key (p,p-1), d = 2p-1, normalized by rho_{p,p-1}(0) (real sqrt2-mult; conj rho = rho)
# LHS: e(1,0)/rho(0) = -i(d-2)
d = 2*p-1
LHS = -i*(d-2)
# parent (p+1,p-1) off-diag: e(0,0)=rho: ratio rho_{p+1,p-1}(0)/rho_{p,p-1}(0):
# K_{d+1}/K_d * (-1) * binom(d+1,p+1)/binom(d,p) = (-i(2d-1)/(d+1))*(-1)*((d+1)/(p+1)) = i(2d-1)/(p+1)
rA0 = i*(2*d-1)/(p+1)          # rho_{p+1,p-1}(0)/rho_{p,p-1}(0)
rA1 = -i*(d-1)*rA0             # e(1,0) ratio
# diagonal parent (p,p): ehat(0,0) = h_hat_p = rho * t10*(4p-3)/p  (derived)
ehat00 = t10*(4*p-3)/p
# ehat(1,0) = (p+1)[t10 e_{p+1,p}(1,0)/(4p-2) + t11 e_{p+1,p}(0,0)/(4p-1)]; rho_{p+1,p}(0)/rho_{p,p-1}(0):
# K_{d+2}/K_d: (-i(2d-1)/(d+1))*(-i(2d+3)/(d+2)) ; binom(d+2,p+1)/binom(d,p) = (d+1)(d+2)/((p+1)(p+2))
# times (-1)^{p+1}/(-1)^p = -1 ; and rho_{p+1,p} has mult 2 like rho_{p,p-1}
rB0 = (-i*(2*d-1)/(d+1))*(-i*(2*d+3)/(d+2)) * (-1) * ((d+1)*(d+2))/((p+1)*(p+2))
rB1 = -i*(d-0)*rB0   # e_{p+1,p}(1,0)/rho_{p,p-1}(0) = -i((d+2)-2)*rho_{p+1,p}(0)/... = -i(d)*rB0
ehat10 = (p+1)*( t10*rB1/(4*p-2) + t11*rB0/(4*p-1) )
# direct DL2 feed: p*[t20*ehat10/(2d-2) + t21*ehat00/(2d-1)]
feed_DL2 = p*( t20*ehat10/(2*d-2) + t21*ehat00/(2*d-1) )
# fold at delta=1: p[conj(g_{4p-5}) + (4p-6) conj(g_{4p-4})], g from ehat with tau1
g_m1 = -( ehat00*t11/(4*p-1) + ehat10*t10/(4*p-2) )   # mu = 4p-5
g_m0 = -( ehat00*t10/(4*p-1) )                        # mu = 4p-4
# conj: t1->t2, and rho real; ehat00 conj = t20*(4p-3)/p; ehat10 conj: substitute t10->t20,t11->t21
ehat10c = (p+1)*( t20*sp.conjugate(rB1)/(4*p-2) + t21*sp.conjugate(rB0)/(4*p-1) )
# rB0, rB1 are Q(i)-rational: conjugate sends i->-i
g_m1c = -( sp.conjugate(ehat00)*t21/(4*p-1) + ehat10c*t20/(4*p-2) )
g_m0c = -( sp.conjugate(ehat00)*t20/(4*p-1) )
feed_fold = p*( g_m1c + (4*p-6)*g_m0c )
# DL1 feed from (p+1,p-1): (p+1)[t10*rA1/(2d-2) + t11*rA0/(2d-1)]
feed_DL1 = (p+1)*( t10*rA1/(2*d-2) + t11*rA0/(2*d-1) )
RHS = feed_DL1 + feed_DL2 + feed_fold
print('simplifying E1-adjacent...')
print(sp.cancel(sp.together(LHS - RHS)))
