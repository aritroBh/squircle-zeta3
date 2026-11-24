
import sympy as sp

sig = sp.Symbol('sig', positive=True)
Z   = sp.Symbol('Z')
tau = sp.Symbol('tau')
sqZ = sp.sqrt(Z)
ORD = 6

y = 1 - sqZ/2*sig + Z/8*sig**2 - Z**2/128*sig**4 + Z**3/1024*sig**6
xs = sp.series(sp.expand(sqZ*sig**-1*y), sig, 0, ORD).removeO()
logy = sp.series(sp.log(1+(y-1)), sig, 0, ORD).removeO()
phi_s = sp.series(sp.expand(sqZ*sig**-1*y - 2*sig**-2*logy), sig, 0, ORD).removeO()
Efac = sp.series(sp.exp(sp.expand(-phi_s + 2*sqZ/sig - Z/2)), sig, 0, ORD).removeO()

x = sp.Symbol('x'); Nsig = sig**-2
def phider_at_xs(j):
    expr = Nsig*(-1)**(j-1)*sp.factorial(j-1)*((xs+Z)**(-j) - xs**(-j))
    return sp.series(expr, sig, 0, ORD).removeO()
phi2 = phider_at_xs(2)
phi_ = {j: phider_at_xs(j) for j in range(3,9)}

def bracket_for_m(mval):
    mm = sp.Rational(2*mval+3, 2)
    # log h ratio series in tau
    lh = mm*(sp.log(1 + tau/(xs+Z))) - sp.log(1 + tau/xs)
    lh = sp.series(lh, tau, 0, 9).removeO()
    lh = sp.series(lh, sig, 0, ORD).removeO()
    Hfac = sp.series(sp.exp(lh), tau, 0, 9).removeO()
    E = sum(-phi_[j]*tau**j/sp.factorial(j) for j in range(3,9))
    expE = sp.series(sp.exp(E), tau, 0, 9).removeO()
    integr = sp.series(sp.expand(expE*Hfac), tau, 0, 9).removeO()
    integr = sp.expand(integr)
    br = 0
    for pw in range(0,9,2):
        c = integr.coeff(tau, pw)
        r = pw//2
        br += sp.expand(c*sp.factorial2(2*r-1)/phi2**r)
    br = sp.series(br, sig, 0, 4).removeO()
    return sp.expand(br)

for mval in [0,1,2]:
    br = bracket_for_m(mval)
    print("m=%d bracket = %s" % (mval, br))
    print()

# ---------- Assembly ----------
# I_s = gam_s*(-i/pi)*e^{-2 sqrt(NZ)} * e^{Z/2} * (pi/2) * sqrt(Z) * N^{-1} * [1 + b1 sig + b2 sig^2 + ...]
# relative weights: m-th local term ~ g_m * Gamma(m+3/2)/Gamma(3/2) * Z^{m/2} sig^m * (series_m)
g1, g2 = sp.symbols('g1 g2')

def series_m(mval):
    mm = sp.Rational(2*mval+3,2)
    # y^{-(m+5/2)} * (y + sig sqZ/2)^{-1/2} * bracket_m * Efac * S_m ; series in sig
    ypow = sp.series(y**(-(sp.Rational(2*mval+5,2))), sig, 0, ORD).removeO()
    Gfac = sp.series((y + sig*sqZ/2)**(-sp.Rational(1,2)), sig, 0, ORD).removeO()
    br = bracket_for_m(mval)
    # Stirling S_m: Gamma(N-m-3/2)/Gamma(N) = N^{-m-3/2}(1 + c1 sig^2 + c2 sig^4)
    c1 = sp.Rational((2*mval+3)*(2*mval+5),8)
    c2 = sp.Rational((2*mval+3)*(2*mval+5)*(2*mval+7)*(6*mval+11),384)
    S = 1 + c1*sig**2 + c2*sig**4
    tot = sp.expand(ypow*Gfac*br*Efac*S)
    tot = sp.series(tot, sig, 0, 4).removeO()
    return sp.expand(tot)

s0 = series_m(0); s1 = series_m(1); s2 = series_m(2)
print("series_m(0) =", s0)
print("series_m(1) =", s1)
print("series_m(2) =", s2)
full = s0 + g1*sp.Rational(3,2)*sqZ*sig*s1 + g2*sp.Rational(15,4)*Z*sig**2*s2
full = sp.expand(sp.series(full, sig, 0, 3).removeO())
b1 = sp.expand(full.coeff(sig,1))
b2 = sp.expand(full.coeff(sig,2))
print("beta1 =", b1)
print("beta2 =", b2)
