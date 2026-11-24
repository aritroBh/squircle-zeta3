# g9_edgelaw5.py -- Proposition 9 (order five): hand route made exact.
#  (A) Psi(u,v) = F^{-2} R on a local edge, to total degree 5 in (uhat, vhat), exact in Q[t, alpha, beta, u0, m]:
#        F = 1 + A uhat + B vhat,  D = alpha uhat - beta vhat,  L F - sigma D = L - alpha sigma0 uhat + beta sigma1 vhat,
#        R = L sqrt((1+Lam)^2 + M),  Lam = (-alpha u0 uhat + beta (u0+1) vhat),  M = m (alpha uhat - beta vhat)^2
#      with u0 = sigma0/L, m = h^2/L^2, A = alpha t, B = beta (1-t).
#      g_ab = [u^a v^b] Psi (uhat = e^u - 1, vhat = e^v - 1).
#  (B) Z5 = 2880 * int sum_{a+b=5, a,b>=1} g_ab,  Z23 = -1440 * int (g32 + g23)   (Beta rules of g9_beta5.py)
#      compared with the formulas of Proposition 9 (xi0=alpha, xi1=beta, L0=-sigma0=-u0 L, L1=L(u0+1)).
#  (C) parallelogram short forms, the regular 2m-gon sums C_m^(5), C_m^(2,3) (symbolic in s, c^2=1-s^2).
import sympy as sp, time
from sympy import Rational as R
t, al, be, u0, m, U, V = sp.symbols('t al be u0 m U V')
GENS = (U, V)
N = 5
def trunc(p):
    p = sp.Poly(p, U, V)
    return sp.Poly.from_dict({k: c for k, c in p.as_dict().items() if sum(k) <= N}, U, V, domain=p.domain) if p.as_dict() else p
def P(x): return sp.Poly(x, U, V, domain=sp.QQ[t, al, be, u0, m])
def mul(a, b): return trunc(a*b)
t0 = time.time()
A = al*t; B = be*(1-t)
w = P(A*U + B*V)
Lam = P(-al*u0*U + be*(u0+1)*V)
M = P(m*(al*U - be*V)**2)
X = P(2*Lam.as_expr()) + mul(Lam, Lam) + M          # (1+Lam)^2 + M - 1
# sqrt(1+X) and (1+w)^(-2) as truncated series
def series(coef, base):
    out = P(1); pw = P(1)
    for k in range(1, N+1):
        pw = mul(pw, base); out = out + pw*coef(k)
    return out
S = series(lambda k: sp.binomial(R(1, 2), k), X)
Finv2 = series(lambda k: sp.binomial(-2, k), w)
Psi = mul(S, Finv2)                                   # Psi / L, as polynomial in uhat, vhat
print('Psi built (%.1f s), %d monomials' % (time.time()-t0, len(Psi.as_dict())))
# convert uhat^i vhat^j -> u^a v^b : [u^a](e^u-1)^i = i! S2(a,i)/a!
from sympy.functions.combinatorial.numbers import stirling
def cc(a, i): return sp.factorial(i)*stirling(a, i)/sp.factorial(a) if i > 0 else (1 if a == 0 else 0)
g = {}
for (i, j), c in Psi.as_dict().items():
    for a in range(i, N+1):
        for b in range(j, N+1-a):
            g[(a, b)] = g.get((a, b), 0) + c*cc(a, i)*cc(b, j)
def integ(expr): return sp.integrate(sp.expand(expr), (t, 0, 1))
g5 = sum(g.get((a, 5-a), 0) for a in range(1, 5))
Z5 = sp.expand(2880*integ(g5))                        # = Z_e^(5) / L
Z23 = sp.expand(-1440*integ(g.get((3, 2), 0) + g.get((2, 3), 0)))
# the law of Proposition 9
xb = -al*u0 + be*(u0+1)                               # xibar = (xi0 L0 + xi1 L1)/L
xp, xt = al+be, al*be
LL = -u0*(u0+1)                                       # L0 L1 / L^2
law5 = 360*xt*m*(20*xb**3 - 32*xb**2 + 15*xb - 2 - 4*(5*xb-4)*xt*LL
                 + m*(5*xt*xp + 8*xp**2 - 28*xt - 15*(xp**2 - 3*xt)*xb))
law23 = 120*xt*m*(12*xb**2 - 12*xb + 2 - 12*(5*xb-4)*xt*LL + m*(15*xt*xp - 3*xp**2 + 24*xt - 45*xt*xb))
d5 = sp.expand(Z5 - law5); d23 = sp.expand(Z23 - law23)
print('Z5/L - law  =', d5, '   Z23/L - law =', d23, '   (%.1f s)' % (time.time()-t0))
assert d5 == 0 and d23 == 0
# sanity: the degree-2 and degree-3 analogues reproduce Z_e of Proposition 7 (zeta(3) part)
g3 = sum(g.get((a, 3-a), 0) for a in (1, 2))
Z3 = sp.expand(12*integ(g3))  # 6 * 2 * int (g21+g12) ; (g30,g03 carry no zeta(3))
law3 = sp.expand(6*al*be*m*(3*xb - 2))
print('order three: Z_e/L - 6 xi0 xi1 (h^2/L^2)(3 xibar - 2) =', sp.expand(Z3 - law3))
# (C) parallelograms: alpha = beta = 1
par5 = sp.expand((Z5).subs({al: 1, be: 1}) - 360*m*(1 + 4*u0*(u0+1) - m))
par23 = sp.expand((Z23).subs({al: 1, be: 1}) - 120*m*(2 + 12*u0*(u0+1) - 3*m))
print('parallelogram short forms, difference:', par5, par23)
# regular 2m-gon: alpha = beta = 2 s^2, L0 = L1 = L/2 (u0 = -1/2), h^2/L^2 = c^2/(4 s^2), h^2/L = c/2, c^2 = 1 - s^2
s, c, mm = sp.symbols('s c mm')
vals = {al: 2*s**2, be: 2*s**2, u0: -R(1, 2), m: (1-s**2)/(4*s**2)}
L_over_hh = None
Z5e = Z5.subs(vals)      # = Z_e^(5)/L ; Z_e^(5) = (Z5e) * L with L = 2 s^2 / c
Ze5 = sp.simplify(Z5e*2*s**2)/c          # times L = 2 s^2/c
Ze23 = sp.simplify(Z23.subs(vals)*2*s**2)/c
C5 = 2880*mm*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)
C23 = -960*mm*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)
d_5 = sp.simplify((2*mm*Ze5 - C5).subs(c, sp.sqrt(1-s**2)))
d_23 = sp.simplify((2*mm*Ze23 - C23).subs(c, sp.sqrt(1-s**2)))
print('regular 2m-gon: 2m Z_e^(5) - C_m^(5) =', d_5, ';  2m Z_e^(2,3) - C_m^(2,3) =', d_23)
e5 = sp.simplify((Ze5 - 1440*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)).subs(c, sp.sqrt(1-s**2)))
e23 = sp.simplify((Ze23 + 480*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)).subs(c, sp.sqrt(1-s**2)))
print('per edge: Z_e^(5) - 1440 s^4 c(...) =', e5, ';  Z_e^(2,3) + 480 c s^4(...) =', e23)
print('done %.1f s' % (time.time()-t0))
