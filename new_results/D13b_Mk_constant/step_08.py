
import sympy as sp

# ---- symbolic saddle expansion, verifying the hand computation ----
# Contribution from z_s:  M_k^(s) = (i s sqrt2 / pi) * ∫_0^inf e^{-u} u^{-1} * J_s(u) du
#   J_s(u) = ∫_0^inf (zeta^{1/2} + c3 zeta^{3/2}) exp(-N log(1+(zs+zeta)/u)) dzeta ,  c3=-1/2, N=k+1
# Expand in eps = N^{-1/2}.  Do it by: log(1+(zs+zeta)/u) = log(1+zs/u) + zeta/(u+zs) - zeta^2/(2(u+zs)^2)+...
# lambda = N/(u+zs). Inner = Gamma(3/2) lam^{-3/2} [1 + (3 c3/2) lam^{-1} + 15/(8N) + O(N^-3/2)]
# Then u-integral of e^{phi(u)} P(u), phi = -u - N log(1+zs/u), P = u^{-1}(u+zs)^{3/2}.
# We verify:  I_s = Gamma(3/2) N^{-3/2} e^{phi(us)} sqrt(2 pi/(-phi''(us))) P(us) [1 + eps*(beta_s) + O(eps^2)]
# with beta_s = (1/2) sqrt(zs) - zs^{3/2}/12   (since 5/4 + 3*c3/2 = 5/4 - 3/4 = 1/2).

N, u = sp.symbols('N, u', positive=True)
eps = sp.symbols('eps')
zs = sp.symbols('zs')

# saddle: u(us+zs) = N zs.  us as series in eps: let us = sqrt(zs)*(1/eps)*U where N = 1/eps^2
# solve exactly: us = (-zs + sqrt(zs^2 + 4 N zs))/2
us_exact = (-zs + sp.sqrt(zs**2 + 4*zs/eps**2))/2
us_ser = sp.series(us_exact, eps, 0, 3).removeO().expand()
print("us =", sp.simplify(us_ser))   # expect sqrt(zs)/eps - zs/2 + zs^{3/2} eps/8

phi = -u - sp.Rational(1,1)/eps**2 * sp.log(1 + zs/u)
phi_us = sp.series(phi.subs(u, us_exact), eps, 0, 3).removeO().expand()
print("phi(us) =", sp.simplify(phi_us))  # expect -2 sqrt(zs)/eps + zs/2 - zs^{3/2} eps/12

phipp = sp.diff(-u - (1/eps**2)*sp.log(1+zs/u), u, 2)
phipp_us = sp.series(phipp.subs(u, us_exact), eps, 0, 2).removeO().expand()
print("phi''(us) =", sp.simplify(phipp_us))

P = (u+zs)**sp.Rational(3,2)/u
P_us = sp.series(P.subs(u, us_exact), eps, 0, 2).removeO().expand()
print("P(us) =", sp.simplify(P_us))

# assemble: prefactor sqrt(2pi/(-phi'')) * P * [1 + (3 c3/2)(us+zs) eps^2 ]  with c3=-1/2
c3 = sp.Rational(-1,2)
bracket = 1 + sp.Rational(3,2)*c3*(us_exact+zs)*eps**2
gauss = sp.sqrt(2*sp.pi/(-phipp.subs(u, us_exact)))
prod = sp.series(sp.expand(gauss*P.subs(u,us_exact)*bracket), eps, 0, 2).removeO().expand()
print("gauss*P*bracket =", sp.simplify(prod))
# expectation: sqrt(pi) zs^{1/4} eps^{-1/2} * (1 + eps*( (5/4 + 3c3/2) sqrt(zs) ))  = ... with c3=-1/2: (1/2) sqrt(zs) eps
