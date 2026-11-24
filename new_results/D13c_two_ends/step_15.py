
import sympy as sp

# Exact saddle s* of the rotated integral: s(1+s) = 2 i N/pi, s* = -1/2 + sqrt(1/4 + 2 i N/pi)
N = sp.symbols('N', positive=True)
i = sp.I
sstar = -sp.Rational(1,2) + sp.sqrt(sp.Rational(1,4) + 2*i*N/sp.pi)
# Large-N series of s*, phase phi(s*) = i pi s*/2 + N log(s*/(1+s*)), phi''(s*), H(s*)
eps = sp.symbols('eps', positive=True)   # eps = 1/sqrt(N)
s_ser = sp.series(sstar.subs(N, 1/eps**2), eps, 0, 4).removeO().expand()
print("s* =", s_ser)
phi = i*sp.pi*sstar/2 + N*sp.log(sstar/(1+sstar))
phi_ser = sp.series(phi.subs(N,1/eps**2), eps, 0, 4).removeO().expand()
print("phi(s*) =", sp.simplify(phi_ser))
# phi''(s) = -i N (1+2s)/(s^2 (1+s)^2)
phi2 = -i*N*(1+2*sstar)/(sstar**2*(1+sstar)**2)
phi2_ser = sp.series(phi2.subs(N,1/eps**2), eps, 0, 3).removeO().expand()
print("phi''(s*) =", sp.simplify(phi2_ser))
# H(s) ~ sqrt(2 pi)/(2 s^{3/2})(1 + a1/s + ...) ; get a1 from gamma ratio: Gamma(s/2)/Gamma(s/2+3/2) = (s/2)^{-3/2}(1 + b1/s + ...)
s = sp.symbols('s', positive=True)
hser = sp.series(sp.gamma(s/2)/sp.gamma(s/2+sp.Rational(3,2)), s, sp.oo, 3)
print("Gamma ratio series:", hser)
# H(s) = sqrt(pi)/4 * ratio ; H(s)/s * sqrt(2 pi/(-phi'')) * e^{phi} * (2/pi) Im -> check leading = -sqrt(pi)/N e^{-th} cos th
Hs = sp.sqrt(sp.pi)/4 * sp.gamma(sstar/2)/sp.gamma(sstar/2+sp.Rational(3,2))
Hs_ser = sp.series(Hs.subs(N,1/eps**2), eps, 0, 3).removeO().expand()
print("H(s*) =", sp.simplify(Hs_ser))
# assemble saddle value SV = e^{phi} sqrt(2 pi/(-phi'')) rot-independent (at s*) times H(s*)/s*
# leading term
SV = sp.exp(phi_ser)*sp.sqrt(2*sp.pi/(-phi2_ser))*Hs_ser/s_ser
SV_ser = sp.series(SV, eps, 0, 3).removeO().expand()
print("SV series (eps=1/sqrt N):", sp.simplify(SV_ser))
